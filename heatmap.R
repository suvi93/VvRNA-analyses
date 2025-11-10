library(pheatmap)
library(data.table)
library(gridExtra)
library(grid)

# Define sample categories
male_samples <- c("MTS53", "MTS54", "MTS55", "MTS56", "MTS57", "MTS58", "MTS59", "MTS60","MLS39", "MLS40", "MLS41", "MLS42", "MLS43", "MLS44", "MLS45", "MHS31", "MHS32", "MHS33", "MHS34", "MHS35", "MHS36", "MHS37", "MHS38")

female_samples <- c("FOS21", "FOS22", "FOS23", "FOS24", "FOS25", "FOS26", "FOS27", "FOS28", "FOS29", "FOS30", "FLS11", "FLS12", "FLS13", "FLS14", "FLS15", "FLS16", "FLS17", "FLS18", "FLS19", "FLS20", "FHS1", "FHS2", "FHS3", "FHS4", "FHS5", "FHS6", "FHS7", "FHS8", "FHS9", "FHS10")

all_samples <- c(male_samples, female_samples)

# Samples to remove
remove_samples <- c("FHS8", "FLS20", "MLS39")  # P24X0

# Remove unwanted samples from the full sample list
all_samples <- setdiff(all_samples, remove_samples)

# Create annotation data frame
sample_sex <- data.frame(Sex = c(rep("Male", length(male_samples)), rep("Female", length(female_samples))))
rownames(sample_sex) <- c(male_samples, female_samples)
sample_sex <- sample_sex[!rownames(sample_sex) %in% remove_samples, , drop = FALSE]

# File list (ordered as: Autosomal, Chr1, ChrX for gonads, head, legs)
files <- c(
  "MF_gonads_onP24XOM_ADEGs_top25TPM.tsv", "MF_head_onP24XOM_ADEGs_top25TPM.tsv", "MF_legs_onP24XOM_ADEGs_top25TPM.tsv",
  "MF_gonads_onP24XOM_Chr1DEGs_top25TPM.tsv",  "MF_head_onP24XOM_Chr1DEGs_top25TPM.tsv",  "MF_legs_onP24XOM_Chr1DEGs_top25TPM.tsv",
  "MF_gonads_onP24XOM_XDEGs_top25TPM.tsv",     "MF_head_onP24XOM_XDEGs_top25TPM.tsv",     "MF_legs_onP24XOM_XDEGs_top25TPM.tsv"
)

# Plot titles
titles <- c("Gonads - Autosomal", "Head - Autosomal", "Legs - Autosomal",
            "Gonads - Chr1", "Head - Chr1", "Legs - Chr1",
            "Gonads - ChrX", "Head - ChrX", "Legs - ChrX")

# Read and prepare matrix: reorder samples to match sex groupings
read_and_prepare <- function(file) {
  df <- fread(file)
  sample_cols <- intersect(all_samples, colnames(df))
  mat <- as.matrix(df[, ..sample_cols])
  mat <- log2(mat + 1)
  rownames(mat) <- df$gene_id
  ordered_cols <- all_samples[all_samples %in% sample_cols]
  mat <- mat[, ordered_cols, drop = FALSE]
  return(mat)
}

# Generate heatmaps
heatmap_list <- list()
for (i in seq_along(files)) {
  mat <- read_and_prepare(files[i])
  annotation <- sample_sex[colnames(mat), , drop = FALSE]

  heatmap <- pheatmap(
    mat,
    silent = TRUE,
    main = titles[i],
    annotation_col = annotation,
    cluster_cols = FALSE,         # Turn off clustering to group by sex explicitly
    show_colnames = FALSE,        # Hide sample names
    fontsize = 7
  )
  heatmap_list[[i]] <- heatmap[[4]]
}

pdf("P24XO_clustered_heatmaps.pdf", width = 12, height = 18)
grid.arrange(grobs = heatmap_list, ncol = 3)
dev.off()

###############################################################################################

library(ggplot2)
library(pheatmap)
library(data.table)
library(gridExtra)
library(grid)

# Define sample categories
male_samples <- c("MTS98", "MTS99", "MTS100", "MTS101", "MTS102", "MTS103", "MTS104", "MLS93", "MLS94", "MLS95", "MLS96", "MLS97", "MHS86", "MHS87", "MHS88", "MHS89", "MHS90", "MHS91", "MHS92")

female_samples <- c("FOS51", "FOS52", "FOS77", "FOS78"," FOS79"," FOS80"," FOS81", "FOS82", "FOS83", "FOS84", "FOS85", "FLS67", "FLS68", "FLS69", "FLS70", "FLS71", "FLS72", "FLS73", "FLS74", "FLS75", "FLS76", "FHS46", "FHS47", "FHS48", "FHS49", "FHS50", "FHS61", "FHS62", "FHS63", "FHS64", "FHS65", "FHS66")

all_samples <- c(male_samples, female_samples)

# Samples to remove
remove_samples <- c("FHS61", "FLS74", "FOS84")  # P24XY

# Remove unwanted samples from the full sample list
all_samples <- setdiff(all_samples, remove_samples)

# Create annotation data frame
sample_sex <- data.frame(Sex = c(rep("Male", length(male_samples)), rep("Female", length(female_samples))))
rownames(sample_sex) <- c(male_samples, female_samples)
sample_sex <- sample_sex[!rownames(sample_sex) %in% remove_samples, , drop = FALSE]

# File list: ordered [Autosomal, Chr1, ChrX] × [Gonads, Head, Legs]
files <- c(
  "MF_gonads_onP24XYF_ADEGs_top25TPM.tsv",      "MF_legs_onP24XYF_ADEGs_top25TPM.tsv",      "MF_head_onP24XYF_ADEGs_top25TPM.tsv",
  "MF_gonads_onP24XYF_XlDEGs_top25TPM.tsv",     "MF_legs_onP24XYF_XlDEGs_top25TPM.tsv",     "MF_head_onP24XYF_XlDEGs_top25TPM.tsv",
  "MF_gonads_onP24XYF_XrYSLRDEGs_top25TPM.tsv", "MF_legs_onP24XYF_XrYSLRDEGs_top25TPM.tsv", "MF_head_onP24XYF_XrYSLRDEGs_top25TPM.tsv",
  "MF_gonads_onP24XYF_XrYPARDEGs_top25TPM.tsv", "MF_legs_onP24XYF_XrYPARDEGs_top25TPM.tsv", "MF_head_onP24XYF_XrYPARDEGs_top25TPM.tsv"
)

titles <- c(
  "Gonads - Autosomal", "Head - Autosomal", "Legs - Autosomal",
  "Gonads - ChrXL",      "Head - ChrXL",      "Legs - ChrXL",
  "Gonads - ChrXR-Y SLR",      "Head - ChrXR-Y SLR",      "Legs - ChrXR-Y SLR",
  "Gonads - ChrXR-Y PAR",      "Head - ChrXR-Y PAR",      "Head - ChrXR-Y PAR" 
)

# Function to check if file is empty
is_file_empty <- function(file) {
  file.info(file)$size == 0
}

# Read and prepare matrix: reorder samples to match sex groupings
read_and_prepare <- function(file) {
  df <- fread(file)
  sample_cols <- intersect(all_samples, colnames(df))
  mat <- as.matrix(df[, ..sample_cols])
  mat <- log2(mat + 1)
  rownames(mat) <- df$gene_id
  ordered_cols <- all_samples[all_samples %in% sample_cols]
  mat <- mat[, ordered_cols, drop = FALSE]
  return(mat)
}

# Generate heatmaps
heatmap_list <- list()
for (i in seq_along(files)) {
 if ((is_file_empty(files[i]))) {
   empty_plot <- ggplot() +
      theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = paste(titles[i], "\n(No data available)"),
               size = 4, hjust = 0.5, vjust = 0.5, fontface = "italic")
    heatmap_list[[i]] <- empty_plot
} else {
  mat <- read_and_prepare(files[i])
  annotation <- sample_sex[colnames(mat), , drop = FALSE]

  heatmap <- pheatmap(
    mat,
    silent = TRUE,
    main = titles[i],
    annotation_col = annotation,
    cluster_cols = FALSE,         # Turn off clustering to group by sex explicitly
    show_colnames = FALSE,        # Hide sample names
    fontsize = 7
  )
  heatmap_list[[i]] <- heatmap[[4]]
}
}

# Save all 9 plots in a 3-column layout (Gonads, Head, Legs per row)
pdf("P24XY_clustered_heatmaps.pdf", width = 18, height = 18)
grid.arrange(grobs = heatmap_list, ncol = 3)
dev.off()
