# ------------------------------------------------------------------
# deseq2.R
#
# Differential gene expression analysis (DESeq2) between males and
# females, run separately per tissue (gonads/head/legs) and per race
# (P24X0/P24XY).
#
# NOTE: this script is a per-tissue template — the count/info file
# names, sample exclusions, and output file names are edited by hand
# for each tissue (gonads/head/legs) and race (P24X0/P24XY) before
# each run, rather than looping over them automatically.
#
# Input:  featureCounts count tables + sample info files from
#         alignment.sh (e.g. MF_gonads_FCforDS.tsv, MF_gonads_info.txt)
# Output: per-tissue/race DEG tables (all genes, male-biased,
#         female-biased) used by heatmap.R and the box_plot_fig*.R
#         scripts
#
# Part of: VvRNA-analyses pipeline (step 2 — differential expression)
# Author:  Suvratha Jayaprasad
# Paper:   Jayaprasad et al. 2026, Genome Biology and Evolution,
#          https://doi.org/10.1093/gbe/evag026
# License: MIT — see LICENSE
# ------------------------------------------------------------------

library(DESeq2)

# ====================================================================
# LOAD COUNT DATA (edit block below for the tissue being processed)
# ====================================================================

count_data = as.matrix(read.csv("MF_gonads_FCforDS.tsv", sep = "\t", row.names = "gene_id"))
count_data <- count_data[, !colnames(count_data) %in% "FOS84"]  # XY
col_data = read.table(file = "MF_gonads_info.txt", header = T, sep = "\t")

count_data = as.matrix(read.csv("MF_head_FCforDS.tsv", sep = "\t", row.names = "gene_id"))
count_data <- count_data[, !colnames(count_data) %in% "FHS8"]   # X0
count_data <- count_data[, !colnames(count_data) %in% "FHS61"]  # XY
col_data = read.table(file = "MF_head_info.txt", header = T, sep = "\t")

count_data = as.matrix(read.csv("MF_legs_FCforDS.tsv", sep = "\t", row.names = "gene_id"))
remove_samples <- c("FLS20", "MLS39")  # X0
count_data <- count_data[, !colnames(count_data) %in% remove_samples]
count_data <- count_data[, !colnames(count_data) %in% "FLS74"]  # XY
col_data = read.table(file = "MF_legs_info.txt", header = T, sep = "\t")

# ====================================================================
# DIFFERENTIAL EXPRESSION
# ====================================================================

dds = DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "FOS")  # changed to FHS for head samples & FLS for leg samples respectively
dds = DESeq(dds)

bmab <- sapply(levels(dds$condition), function(lvl) rowMeans(counts(dds, normalized = TRUE)[, dds$condition == lvl]))  # to get base expression of male and female samples
res <- results(dds, alpha = 0.05, lfcThreshold = 2)  # padj and log2FC cutoffs
z <- cbind(gene_id = rownames(res), bmab, res)
z$padj[!is.na(z$pvalue) & z$pvalue == 1] <- 1  # to ensure genes are not treated as false-positives

z_filtered_male <- subset(z, padj < 0.05 & log2FoldChange > 2)    # male-biased genes
z_filtered_female <- subset(z, padj < 0.05 & log2FoldChange < -2)  # female-biased genes

# ====================================================================
# WRITE OUTPUT (edit file names below for the tissue/race being processed)
# ====================================================================

# for P24X0 race, changed gonads to head / legs respectively
write.table(z, file = "MF_gonads_onP24XOM_DEGs.tsv", sep = '\t', quote = F, row.names = F)
write.table(z_filtered_male, file = "MF_gonads_onP24XOM_mbDEGs.tsv", sep = '\t', quote = F, row.names = F)
write.table(z_filtered_female, file = "MF_gonads_onP24XOM_fbDEGs.tsv", sep = '\t', quote = F, row.names = F)

# for P24XY race, changed gonads to head / legs respectively
write.table(z, file = "MF_gonads_onP24XYF_DEGs.tsv", sep = '\t', quote = F, row.names = F)
write.table(z_filtered_male, file = "MF_gonads_onP24XYF_mbDEGs.tsv", sep = '\t', quote = F, row.names = F)
write.table(z_filtered_female, file = "MF_gonads_onP24XYF_fbDEGs.tsv", sep = '\t', quote = F, row.names = F)
