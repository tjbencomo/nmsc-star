# Description: Perform differential expression tests on NMSC data quantified using STAR-RSEM

library(DESeq2)
library(sva)
library(tximport)
library(IHW)
library(apeglm)
library(readr)
library(dplyr)
library(stringr)
library(BiocParallel)
register(MulticoreParam(8))

# Setup
rsemDir <- "rsem"
# This dir will be named scc-meta-analysis if you are running scripts
# after reading the paper
dataDir <- file.path("~", "nmsc-rna-seq", "data")
metadataFp <- file.path(dataDir, "metadata_final_cohort.csv")
tx2gFp <- file.path(dataDir, "gencode.v38.tx2gene.csv.gz")
deseqDir <- "deseq"
if (!dir.exists(deseqDir)) {
    dir.create(deseqDir)
}

# Load metadata
print("Loading metadata")
metadata <- read_csv(metadataFp)
tx2g <- read_csv(tx2gFp)
condition_order <- c("NS", "AK", "AK_IEC", "AK_IEC_SCC", "IEC", "KA", "SCC")

# metadata <- metadata %>%
#     filter(condition != "BCC") %>%
#     mutate(condition = factor(condition, levels = condition_order)) %>%
#     mutate(study = factor(study_name)) %>%
#     as.data.frame()


# print(table(metadata$condition))
# print(table(metadata$study))
# print(dim(metadata))

# # Load RSEM quant files
# print("Loading RSEM files")
# files <- file.path(rsemDir, paste0(metadata$Sample, ".isoforms.results"))
# names(files) <- metadata$Sample
# txi.rsem <- tximport(files, type = "rsem", txIn = TRUE, txOut = FALSE, tx2gene = tx2g %>% select(-gene_symbol))


# # Check that samples are properly matched with read data
# stopifnot(all(colnames(txi.rsem$counts) == metadata$Sample))
# rownames(metadata) <- metadata$Sample

# print("Creating DESeq2 object")
# dds <- DESeqDataSetFromTximport(txi.rsem , metadata, ~ study + condition)
# # dds <- DESeqDataSetFromTximport(txi.rsem , metadata, ~ condition)
# print("Running DESeq2")
# dds <- DESeq(dds, parallel = T)

# print(dds)
# saveRDS(dds, "deseq/tmp.deseq.rds")
dds <- readRDS("deseq/tmp.deseq.rds")

# Add gene symbol info to dds object
genedf <- tx2g %>% select(-transcript_id) %>% distinct(gene_id, .keep_all = T)
genedf <- genedf[match(rownames(dds), genedf$gene_id), ]
genedf$ensgene <- str_split(genedf$gene_id, "\\.", simplify = T)[, 1]
stopifnot(all(rownames(dds) == genedf$gene_id))
rowData(dds)$ensgene <- genedf$ensgene
rowData(dds)$symbol <- genedf$gene_symbol

# Compute DEGs
contrast_list <- list(
    'NS_vs_AK' = c('NS', 'AK'),
    'AK_vs_SCC' = c('AK', 'SCC'),
    'NS_vs_SCC' = c('NS', 'SCC'),
    'NS_vs_IEC' = c('NS', 'IEC'),
    'NS_vs_KA' = c('NS', 'KA'),
    'IEC_vs_KA' = c('IEC', 'KA'),
    'AK_vs_IEC' = c('AK', 'IEC'),
    'AK_vs_IEC' = c('AK', 'KA'),
    'IEC_vs_SCC' = c('IEC', 'SCC'),
    'KA_vs_SCC' = c('KA', 'SCC')
)

alphaThresh <- .05
lfcThresh <- .5
print("Computing contrasts")
for (i in 1:length(contrast_list)) {
    contrast_name <- names(contrast_list)[i]
    print(paste("Computing DEGs for", contrast_name))
    numerator_name <- contrast_list[[i]][2]
    denominator_name <- contrast_list[[i]][1]
    direction_col_name <- paste0('up_in_', numerator_name)
    print(direction_col_name)
    ## Get and save thresholded test results
    print("Calculating lfcThreshold=.5 DEGs")
    mle_res_thresh <- results(dds, filterFun = ihw, alpha = alphaThresh, 
                       contrast = c("condition", numerator_name, denominator_name),
                       parallel = T, lfcThreshold = lfcThresh)
    mle_res_thresh$gene_symbol <- rowData(dds)$symbol
    mle_res_thresh_df <- mle_res_thresh %>% as.data.frame() %>% tibble::rownames_to_column('gene_id') %>%
               mutate(!!direction_col_name := log2FoldChange > 0) 
    fn <- paste0(numerator_name, '_vs_', denominator_name, "_lfcThreshold", '.csv') 
    fp <- file.path(deseqDir, fn)
    print(fp)
    write_csv(mle_res_thresh_df, fp)
   
    ## Get and save non-thresholded results
    print("Calculating non-thresholded DEGs")
    mle_res <-results(dds, filterFun = ihw, alpha = alphaThresh, 
                      contrast = c("condition", numerator_name, denominator_name),
                      parallel = T )
    mle_res$gene_symbol <- rowData(dds)$symbol 
    mle_res_df <- mle_res %>% as.data.frame() %>% tibble::rownames_to_column('gene_id') %>%
               mutate(!!direction_col_name := log2FoldChange > 0)
    fn2 <- paste0(numerator_name, '_vs_', denominator_name, '.csv') 
    fp2 <- file.path(deseqDir, fn2)
    print(fp2)
    write_csv(mle_res_df, fp2)

    ## Get and save MAP results if available
    coef_name <- paste0("condition_", numerator_name, "_vs_", denominator_name)
    if (coef_name %in% resultsNames(dds)) {
        print("Calculating MAP DEGs")
        map_res <- lfcShrink(dds, coef=coef_name, 
                             type="apeglm", parallel = T)
        map_res$gene_symbol <- rowData(dds)$symbol
        map_res_df <- map_res %>% as.data.frame() %>% tibble::rownames_to_column('gene_id') %>%
                   mutate(!!direction_col_name := log2FoldChange > 0)
        fn3 <- paste0(numerator_name, '_vs_', denominator_name, "_MAP", '.csv') 
        fp3 <- file.path(deseqDir, fn3)
        print(fp3)
        write_csv(map_res_df, fp3)
    } else {
        print(paste("Coefficient", coef_name, "not found in resultsNames(dds)"))
        print("Skipping MAP estimation")
    }
}

contrast_list_map <- contrast_list[c("NS_vs_AK", "AK_vs_SCC", "NS_vs_SCC")]
print("Calculating MAP estimates")

# VST normalization without batch correction
vsd <- vst(dds, blind = F)
vsdBlind = vst(dds, blind = T)

# Batch correction
# See batch_correction.R in nmsc-rna-seq repo

# Save results
print("Saving results")
saveRDS(dds, file.path(deseqDir, "deseq_obj.rds"))
saveRDS(vsd, file.path(deseqDir, "vst_normalized_counts.rds"))
saveRDS(vsdBlind, file.path(deseqDir, "vst_normalized_counts_blind.rds"))

