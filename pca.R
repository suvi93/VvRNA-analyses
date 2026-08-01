# ------------------------------------------------------------------
# pca.R
#
# Principal component analysis of rlog-transformed expression data,
# for P24X0 and P24XY, with and without outlier samples removed.
#
# Input:  per-race, all-tissue featureCounts tables + sample info
#         files (e.g. MF_alltissues_FC_P24X0.tsv)
# Output: 4-panel PCA figure (All_Tissues_PCA.pdf)
#
# Part of: VvRNA-analyses pipeline (exploratory analysis / QC)
# Author:  Suvratha Jayaprasad
# Paper:   Jayaprasad et al. 2026, Genome Biology and Evolution,
#          https://doi.org/10.1093/gbe/evag026
# License: MIT — see LICENSE
# ------------------------------------------------------------------

library(DESeq2)
library(ggplot2)
library(ggpubr)
library(patchwork)

# ====================================================================
# PCA PLOT FUNCTION
# ====================================================================

# Function to generate PCA plot for a dataset
generate_pca_plot <- function(count_file, info_file, title, outlier_samples = NULL) {

  # Step 1: Read data
  count_data <- read.csv(count_file, sep = "\t", row.names = "gene_id")
  col_data <- read.table(info_file, header = TRUE, sep = "\t", row.names = 1)

  # Step 2: Ensure columns and metadata align
  stopifnot(all(colnames(count_data) %in% rownames(col_data)))
  col_data <- col_data[colnames(count_data), , drop = FALSE]

  # Step 3: Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = col_data,
                                design = ~ condition)

  # Step 4: vst transformation
  rld <- rlogTransformation(dds)

  # Identify and optionally remove outliers
  if (!is.null(outlier_samples)) {
    dds <- dds[, !(colnames(dds) %in% outlier_samples)]
    rld <- rlogTransformation(dds)
  }

  # Step 5: Generate PCA plot
  pca_plot <- plotPCA(rld, intgroup = "condition") +
    ggtitle(title) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11)
    ) +
    scale_color_manual(
      name = "Tissues",
      values = c(
        "MHS" = "#1f77b4",  # Male head
        "FHS" = "#ff7f0e",  # Female head
        "MLS" = "#2ca02c",  # Male legs
        "FLS" = "#d62728",  # Female legs
        "MTS" = "#9467bd",  # Male gonads
        "FOS" = "#8c564b"   # Female gonads
      ),
      labels = c(
        "MHS" = "Male Head", "FHS" = "Female Head",
        "MLS" = "Male Legs", "FLS" = "Female Legs",
        "MTS" = "Male Testes", "FOS" = "Female Ovaries"
      )
    )

  return(pca_plot)
}

# ====================================================================
# GENERATE ALL PANELS
# ====================================================================

# PANEL A — P24X0 (all samples)
p1 <- generate_pca_plot(
  count_file = "MF_alltissues_FC_P24X0.tsv",
  info_file  = "MF_alltissues_info_P24X0.txt",
  title = "(a) P24X0 - All Samples"
)

# PANEL B — P24X0 (outliers removed)
p2 <- generate_pca_plot(
  count_file = "MF_alltissues_FC_P24X0.tsv",
  info_file  = "MF_alltissues_info_P24X0.txt",
  outlier_samples = c("FHS8", "FLS20", "MLS39"),
  title = "(b) P24X0 - Outliers Removed"
)

# PANEL C — P24XY (all samples)
p3 <- generate_pca_plot(
  count_file = "MF_alltissues_FC_P24XY.tsv",
  info_file  = "MF_alltissues_info_P24XY.txt",
  title = "(c) P24XY - All Samples"
)

# PANEL D — P24XY (outliers removed)
p4 <- generate_pca_plot(
  count_file = "MF_alltissues_FC_P24XY.tsv",
  info_file  = "MF_alltissues_info_P24XY.txt",
  outlier_samples = c("FHS61", "FLS74", "FOS84"),
  title = "(d) P24XY - Outliers Removed"
)

# ====================================================================
# COMBINE ALL PANELS
# ====================================================================
combined_plot <- (p1 | p2) / (p3 | p4)

# ====================================================================
# SAVE ALL PANELS
# ====================================================================
ggsave("All_Tissues_PCA.pdf", combined_plot, width = 12, height = 10)
