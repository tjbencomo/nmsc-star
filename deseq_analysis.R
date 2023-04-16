# Description: Perform differential expression tests on NMSC data quantified using STAR-RSEM

library(DESeq2)
library(tximport)
library(IHW)
library(readr)
library(dplyr)
library(stringr)
library(BiocParallel)
register(MulticoreParam(4))

# Setup
rsemDir <- "rsem"
dataDir <- file.path("~", "nmsc-rna-seq", "data")
metadataFp <- file.path(dataDir, "Final_Cohort.csv")
tx2gFp <- file.path(dataDir, "gencode.v38.tx2gene.tsv.gz")
deseqDir <- "deseq"
if (!dir.exists(deseqDir)) {
    dir.create(deseqDir)
}

# Load metadata
print("Loading metadata")
metadata <- read_csv(metadataFp)
tx2g <- read_tsv(tx2gFp) %>% select(-hgnc_symbol)
metadata <- metadata %>%
    filter(condition != "BCC") %>%
    mutate(condition = factor(condition, levels = c('NS', 'IEC', 'Mixed', 'BCC', 'AK', 'SCC'))) %>%
    mutate(study = factor(study_accession)) %>%
    as.data.frame()

# Load RSEM quant files
print("Loading RSEM files")
files <- file.path(rsemDir, paste0(metadata$Sample, ".isoforms.results"))
names(files) <- metadata$Sample
txi.rsem <- tximport(files, type = "rsem", txIn = TRUE, txOut = FALSE, tx2gene = tx2g)


# Check that samples are properly matched with read data
stopifnot(all(colnames(txi.rsem$counts) == metadata$Sample))
rownames(metadata) <- metadata$Sample

print("Creating DESeq2 object")
dds <- DESeqDataSetFromTximport(txi.rsem , metadata, ~ study + condition)
print("Running DESeq2")
dds <- DESeq(dds, parallel = T)

print(dds)

print("Computing contrasts")
resSCC_NS <- results(dds, filterFun = ihw, alpha = .01, contrast = c("condition", "SCC", "NS"), parallel = T)
resSCC_AK <- results(dds, filterFun = ihw, alpha = .01, contrast = c("condition", "SCC", "AK"), parallel = T)
resAK_NS <- results(dds, filterFun = ihw, alpha = .01, contrast = c("condition", "AK", "NS"), parallel = T)
resIEC_NS <- results(dds, filterFun = ihw, alpha = .01, contrast = c("condition", "IEC", "NS"), parallel = T)
resMixed_NS <- results(dds, filterFun = ihw, alpha = .01, contrast = c("condition", "Mixed", "NS"), parallel = T)
resSCC_IEC <- results(dds, filterFun = ihw, alpha = .01, contrast = c("condition", "SCC", "IEC"), parallel = T)

print("Saving results")
saveRDS(dds, file.path(deseqDir, "deseq_obj.rds"))
write_csv(data.frame(resSCC_NS), file.path(deseqDir, "SCC_vs_NS.csv"))
write_csv(data.frame(resSCC_AK), file.path(deseqDir, "SCC_vs_AK.csv"))
write_csv(data.frame(resAK_NS), file.path(deseqDir, "AK_vs_NS.csv"))
write_csv(data.frame(resIEC_NS), file.path(deseqDir, "IEC_vs_NS.csv"))
write_csv(data.frame(resMixed_NS), file.path(deseqDir, "Mixed_vs_NS.csv"))
write_csv(data.frame(resSCC_IEC), file.path(deseqDir, "SCC_vs_IEC.csv"))

# saveRDS(resSCC_NS, file.path(deseqDir, "SCC_vs_NS.rds"))
# saveRDS(resSCC_AK, file.path(deseqDir, "SCC_vs_AK.rds"))
# saveRDS(resAK_NS, file.path(deseqDir, "AK_vs_NS.rds"))
# saveRDS(resIEC_NS, file.path(deseqDir, "IEC_vs_NS.rds"))
# saveRDS(resMixed_NS, file.path(deseqDir, "Mixed_vs_NS.rds"))
# saveRDS(resSCC_IEC, file.path(deseqDir, "SCC_vs_IEC.rds"))


