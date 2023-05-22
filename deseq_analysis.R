# Description: Perform differential expression tests on NMSC data quantified using STAR-RSEM

library(DESeq2)
library(tximport)
library(IHW)
library(apeglm)
library(readr)
library(dplyr)
library(stringr)
library(BiocParallel)
register(MulticoreParam(4))

# Setup
rsemDir <- "rsem"
dataDir <- file.path("~", "nmsc-rna-seq", "data")
metadataFp <- file.path(dataDir, "Final_Cohort.csv")
tx2gFp <- file.path(dataDir, "gencode.v38.tx2gene.csv.gz")
deseqDir <- "deseq"
if (!dir.exists(deseqDir)) {
    dir.create(deseqDir)
}

# Load metadata
print("Loading metadata")
metadata <- read_csv(metadataFp)
tx2g <- read_csv(tx2gFp)
metadata <- metadata %>%
    filter(condition != "BCC") %>%
    mutate(condition = factor(condition, levels = c('NS', 'IEC', 'Mixed', 'BCC', 'AK', 'SCC'))) %>%
    mutate(study = factor(study_accession)) %>%
    as.data.frame()

# Load RSEM quant files
print("Loading RSEM files")
files <- file.path(rsemDir, paste0(metadata$Sample, ".isoforms.results"))
names(files) <- metadata$Sample
txi.rsem <- tximport(files, type = "rsem", txIn = TRUE, txOut = FALSE, tx2gene = tx2g %>% select(-gene_symbol))


# Check that samples are properly matched with read data
stopifnot(all(colnames(txi.rsem$counts) == metadata$Sample))
rownames(metadata) <- metadata$Sample

print("Creating DESeq2 object")
dds <- DESeqDataSetFromTximport(txi.rsem , metadata, ~ study + condition)
print("Running DESeq2")
dds <- DESeq(dds, parallel = T)

print(dds)
# saveRDS(dds, "deseq/tmp.deseq.rds")

# Add gene symbol info to dds object
genedf <- tx2g %>% select(-transcript_id) %>% distinct(gene_id, .keep_all = T)
genedf <- genedf[match(rownames(dds), genedf$gene_id), ]
genedf$ensgene <- str_split(genedf$gene_id, "\\.", simplify = T)[, 1]
stopifnot(all(rownames(dds) == genedf$gene_id))
rowData(dds)$ensgene <- genedf$ensgene
rowData(dds)$symbol <- genedf$gene_symbol

# Compute DEGs
print("Computing contrasts")
alphaThresh <- .01
resSCC_NS <- results(dds, filterFun = ihw, alpha = alphaThresh, contrast = c("condition", "SCC", "NS"), parallel = T) %>%
    data.frame() %>%
    tibble::rownames_to_column('gene_id') %>%
    inner_join(genedf) %>%
    mutate(up_in_scc = log2FoldChange > 0)
resSCC_AK <- results(dds, filterFun = ihw, alpha = alphaThresh, contrast = c("condition", "SCC", "AK"), parallel = T) %>%
    data.frame() %>%
    tibble::rownames_to_column('gene_id') %>%
    inner_join(genedf) %>%
    mutate(up_in_scc = log2FoldChange > 0)
resAK_NS <- results(dds, filterFun = ihw, alpha = alphaThresh, contrast = c("condition", "AK", "NS"), parallel = T) %>%
    data.frame() %>%
    tibble::rownames_to_column('gene_id') %>%
    inner_join(genedf) %>%
    mutate(up_in_ak = log2FoldChange > 0)
resIEC_NS <- results(dds, filterFun = ihw, alpha = alphaThresh, contrast = c("condition", "IEC", "NS"), parallel = T) %>%
    data.frame() %>%
    tibble::rownames_to_column('gene_id') %>%
    inner_join(genedf) %>%
    mutate(up_in_iec = log2FoldChange > 0)
resMixed_NS <- results(dds, filterFun = ihw, alpha = alphaThresh, contrast = c("condition", "Mixed", "NS"), parallel = T) %>%
    data.frame() %>%
    tibble::rownames_to_column('gene_id') %>%
    inner_join(genedf) %>%
    mutate(up_in_mixed = log2FoldChange > 0)
resSCC_IEC <- results(dds, filterFun = ihw, alpha = alphaThresh, contrast = c("condition", "SCC", "IEC"), parallel = T) %>%
    data.frame() %>%
    tibble::rownames_to_column('gene_id') %>%
    inner_join(genedf) %>%
    mutate(up_in_scc = log2FoldChange > 0)
resSCC_Mixed <- results(dds, filterFun = ihw, alpha = alphaThresh, contrast = c("condition", "SCC", "Mixed"), parallel = T) %>%
    data.frame() %>%
    tibble::rownames_to_column('gene_id') %>%
    inner_join(genedf) %>%
    mutate(up_in_scc = log2FoldChange > 0)
resIEC_Mixed <- results(dds, filterFun = ihw, alpha = alphaThresh, contrast = c("condition", "IEC", "Mixed"), parallel = T) %>%
    data.frame() %>%
    tibble::rownames_to_column('gene_id') %>%
    inner_join(genedf) %>%
    mutate(up_in_iec = log2FoldChange > 0)

resIEC_AK <- results(dds, filterFun = ihw, alpha = alphaThresh, contrast = c("condition", "IEC", "AK"), parallel = T) %>%
    data.frame() %>%
    tibble::rownames_to_column('gene_id') %>%
    inner_join(genedf) %>%
    mutate(up_in_iec = log2FoldChange > 0)


# MAP estimates
resMapSCC_NS <- lfcShrink(dds, coef = "condition_SCC_vs_NS", type = "apeglm") %>%
    data.frame() %>%
    tibble::rownames_to_column('gene_id') %>%
    inner_join(genedf) %>%
    mutate(up_in_scc = log2FoldChange > 0)
# coef doesn't exist - this is the weird contrast issue with apeglm vs the other methods?
# resMapSCC_AK <- lfcShrink(dds, coef = "condition_SCC_vs_AK", type = "apeglm") %>%
#     data.frame() %>%
#     tibble::rownames_to_column('gene_id') %>%
#     inner_join(genedf) %>%
#     mutate(up_in_scc = log2FoldChange > 0)
resMapAK_NS <- lfcShrink(dds, coef = "condition_AK_vs_NS", type = "apeglm") %>%
    data.frame() %>%
    tibble::rownames_to_column('gene_id') %>%
    inner_join(genedf) %>%
    mutate(up_in_ak = log2FoldChange > 0)



# VST normalization
vsd <- vst(dds, blind = F)
vsdBlind = vst(dds, blind = T)

# Batch correction
# See batch_correction.R in nmsc-rna-seq repo

# Save results
print("Saving results")
saveRDS(dds, file.path(deseqDir, "deseq_obj.rds"))
saveRDS(vsd, file.path(deseqDir, "vst_normalized_counts.rds"))
saveRDS(vsdBlind, file.path(deseqDir, "vst_normalized_counts_blind.rds"))
write_csv(resSCC_NS, file.path(deseqDir, "SCC_vs_NS.csv"))
write_csv(resSCC_AK, file.path(deseqDir, "SCC_vs_AK.csv"))
write_csv(resAK_NS, file.path(deseqDir, "AK_vs_NS.csv"))
write_csv(resIEC_NS, file.path(deseqDir, "IEC_vs_NS.csv"))
write_csv(resMixed_NS, file.path(deseqDir, "Mixed_vs_NS.csv"))
write_csv(resSCC_IEC, file.path(deseqDir, "SCC_vs_IEC.csv"))
write_csv(resSCC_Mixed, file.path(deseqDir, "SCC_vs_IEC.csv"))
write_csv(resIEC_Mixed, file.path(deseqDir, "SCC_vs_IEC.csv"))
write_csv(resIEC_AK, file.path(deseqDir, "IEC_AK.csv"))
write_csv(resMapSCC_NS, file.path(deseqDir, "MAP_SCC_vs_NS.csv"))
# write_csv(resMapSCC_AK, file.path(deseqDir, "MAP_SCC_vs_AK.csv"))
write_csv(resMapAK_NS, file.path(deseqDir, "MAP_AK_vs_NS.csv"))

