library(readr)
library(dplyr)
library(ggplot2)
library(ggstatsplot)
library(RColorBrewer)
library(patchwork)

setwd("/Users/octaviopalacios/Desktop/DC_paper/Figures/Figure_1")

# Function to load and filter data
load_and_filter_data <- function(file, female_col, male_col, chr_label) {
  df <- read.delim(file, header = TRUE)
  df <- df %>% filter(.data[[female_col]] >= 1, .data[[male_col]] >= 1)
  df$gene_id <- chr_label
  return(df)
}

# Function to prepare dataset for plotting
prepare_data <- function(df, column_name, new_col_name) {
  df <- df %>% select(gene_id, all_of(column_name))
  colnames(df)[2] <- new_col_name
  return(df)
}

# Function to generate boxplot
plot_boxplot <- function(data, value_col, title, palette_used = NULL) {
  data <- data %>% mutate(log2_value = log2(.data[[value_col]]))
  plot <- ggbetweenstats(
    data = data, x = gene_id, y = log2_value,
    xlab = "", ylab = "log2(TPM)",
    type = "np", plot.type = "box",
    p.adjust.method = "BH", conf.level = 0.99,
    pairwise.display = "significant", title = title,
    package = "ggsci", palette = "nrc_npg",
    boxplot.args = list(width = 0.5, alpha = 0.2, na.rm = TRUE),
    violin.args = list(width = 0.7, alpha = 0.2, na.rm = TRUE),
    centrality.label.args = list(size = 2, nudge_x = 0.4, segment.linetype = 4,
                                 min.segment.length = 0)
  )
  if(!is.null(palette_used)){
    plot <- plot + scale_color_manual(name = "Chromosome", values = palette_used)
  }
  plot + theme(legend.position = 'none', legend.text = element_text(size = 10))
}

# Load and process datasets for head and legs
body_parts <- list(
  "head" = list(
    "autosomes" = "MF_head_autosomes_TPMmean.tsv",
    "chr1" = "MF_head_chr1_TPMmean.tsv",
    "chrx" = "MF_head_chrx_TPMmean.tsv"
  ),
  "legs" = list(
    "autosomes" = "MF_legs_autosomes_TPMmean.tsv",
    "chr1" = "MF_legs_chr1_TPMmean.tsv",
    "chrx" = "MF_legs_chrx_TPMmean.tsv"
  )
)

# Initialize lists to store combined data
Male_combined <- list()
Female_combined <- list()

for (part in names(body_parts)) {
  datasets <- body_parts[[part]]
  
  Male_combined[[part]] <- bind_rows(
    prepare_data(load_and_filter_data(datasets$autosomes, paste0("F", toupper(substr(part, 1, 1)), "S_autosomes_mean"), paste0("M", toupper(substr(part, 1, 1)), "S_autosomes_mean"), "Autosomes"), paste0("M", toupper(substr(part, 1, 1)), "S_autosomes_mean"), "Mean_TPM"),
    prepare_data(load_and_filter_data(datasets$chr1, paste0("F", toupper(substr(part, 1, 1)), "S_chr1_mean"), paste0("M", toupper(substr(part, 1, 1)), "S_chr1_mean"), "chr1"), paste0("M", toupper(substr(part, 1, 1)), "S_chr1_mean"), "Mean_TPM"),
    prepare_data(load_and_filter_data(datasets$chrx, paste0("F", toupper(substr(part, 1, 1)), "S_chrx_mean"), paste0("M", toupper(substr(part, 1, 1)), "S_chrx_mean"), "chrX"), paste0("M", toupper(substr(part, 1, 1)), "S_chrx_mean"), "Mean_TPM")
  )
  
  Female_combined[[part]] <- bind_rows(
    prepare_data(load_and_filter_data(datasets$autosomes, paste0("F", toupper(substr(part, 1, 1)), "S_autosomes_mean"), paste0("M", toupper(substr(part, 1, 1)), "S_autosomes_mean"), "Autosomes"), paste0("F", toupper(substr(part, 1, 1)), "S_autosomes_mean"), "Mean_TPM"),
    prepare_data(load_and_filter_data(datasets$chr1, paste0("F", toupper(substr(part, 1, 1)), "S_chr1_mean"), paste0("M", toupper(substr(part, 1, 1)), "S_chr1_mean"), "chr1"), paste0("F", toupper(substr(part, 1, 1)), "S_chr1_mean"), "Mean_TPM"),
    prepare_data(load_and_filter_data(datasets$chrx, paste0("F", toupper(substr(part, 1, 1)), "S_chrx_mean"), paste0("M", toupper(substr(part, 1, 1)), "S_chrx_mean"), "chrX"), paste0("F", toupper(substr(part, 1, 1)), "S_chrx_mean"), "Mean_TPM")
  )
}

# Generate color palette
colourCount <- 14
palette1 <- rev(colorRampPalette(brewer.pal(12, "Paired"))(colourCount))

# Generate plots
p1 <- plot_boxplot(Male_combined$head, "Mean_TPM", "Male Head", palette1)
p2 <- plot_boxplot(Female_combined$head, "Mean_TPM", "Female Head", palette1)
p3 <- plot_boxplot(Male_combined$legs, "Mean_TPM", "Male Legs", palette1)
p4 <- plot_boxplot(Female_combined$legs, "Mean_TPM", "Female Legs", palette1)

# Modify the jitter layer inside the ggplot object
# Identify which layer uses position_jitter and change its width
p1$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)
p2$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)
p3$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)
p4$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)

print(p1)
print(p2)
print(p3)
print(p4)

# Compute Male/Female ratio separately for heads and legs
MF_ratio_head <- Male_combined$head %>% mutate(Ratio = Mean_TPM / Female_combined$head$Mean_TPM) %>% mutate(log2_Ratio = log2(Ratio))
MF_ratio_legs <- Male_combined$legs %>% mutate(Ratio = Mean_TPM / Female_combined$legs$Mean_TPM) %>% mutate(log2_Ratio = log2(Ratio))

# Plot Male/Female ratio separately for heads and legs
p5 <- ggbetweenstats(
  data = MF_ratio_head, x = gene_id, y = log2_Ratio,
  xlab = "", ylab = "log2(M/F TPM)",
  type = "np", plot.type = "box",
  p.adjust.method = "BH", conf.level = 0.99,
  pairwise.display = "significant", title = "Male to Female Ratio - Head",
  package = "ggsci", palette = "nrc_npg",
  boxplot.args = list(width = 0.5, alpha = 0.2, na.rm = TRUE),
  violin.args = list(width = 0.7, alpha = 0.2, na.rm = TRUE),
  centrality.label.args = list(size = 2, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)
) + scale_color_manual(name = "Chromosome", values = palette1) + theme(legend.position = 'none')

p6 <- ggbetweenstats(
  data = MF_ratio_legs, x = gene_id, y = log2_Ratio,
  xlab = "", ylab = "log2(M/F TPM)",
  type = "np", plot.type = "box",
  p.adjust.method = "BH", conf.level = 0.99,
  pairwise.display = "significant", title = "Male to Female Ratio - Leg",
  package = "ggsci", palette = "nrc_npg",
  boxplot.args = list(width = 0.5, alpha = 0.2, na.rm = TRUE),
  violin.args = list(width = 0.7, alpha = 0.2, na.rm = TRUE),
  centrality.label.args = list(size = 2, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)
) + scale_color_manual(name = "Chromosome", values = palette1) + theme(legend.position = 'none')

# Modify the jitter layer inside the ggplot object
# Identify which layer uses position_jitter and change its width
p5$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)
p6$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)

print(p5)
print(p6)

## Foor the gonads
# Function to load and filter data
load_and_filter_data <- function(file, female_col, male_col, chr_label) {
  df <- read.delim(file, header = TRUE, check.names = TRUE)  # Ensure correct column names
  
  # Debug: Print column names to verify
  print(paste("Loading file:", file))
  print("Column names detected:")
  print(colnames(df))
  
  # Check if required columns exist before filtering
  if (!(female_col %in% colnames(df)) | !(male_col %in% colnames(df))) {
    stop(paste("Error: Column", female_col, "or", male_col, "not found in", file))
  }
  
  # Apply filtering only if columns exist
  df <- df %>% filter(.data[[female_col]] >= 1, .data[[male_col]] >= 1)
  df$gene_id <- chr_label  # Label the chromosome group
  return(df)
}

# Function to prepare dataset for plotting
prepare_data <- function(df, column_name, new_col_name) {
  df <- df %>% select(gene_id, all_of(column_name))
  colnames(df)[2] <- new_col_name
  return(df)
}

# Load and process gonad datasets
gonad_files <- list(
  "autosomes" = "MF_gonads_autosomes_TPMmean.tsv",
  "chr1" = "MF_gonads_chr1_TPMmean.tsv",
  "chrx" = "MF_gonads_chrx_TPMmean.tsv"
)

Male_combined <- bind_rows(
  prepare_data(load_and_filter_data(gonad_files$autosomes, "FOS_autosomes_mean", "MTS_autosomes_mean", "Autosomes"), "MTS_autosomes_mean", "Mean_TPM"),
  prepare_data(load_and_filter_data(gonad_files$chr1, "FOS_chr1_mean", "MTS_chr1_mean", "chr1"), "MTS_chr1_mean", "Mean_TPM"),
  prepare_data(load_and_filter_data(gonad_files$chrx, "FOS_chrx_mean", "MTS_chrx_mean", "chrX"), "MTS_chrx_mean", "Mean_TPM")
)

Female_combined <- bind_rows(
  prepare_data(load_and_filter_data(gonad_files$autosomes, "FOS_autosomes_mean", "MTS_autosomes_mean", "Autosomes"), "FOS_autosomes_mean", "Mean_TPM"),
  prepare_data(load_and_filter_data(gonad_files$chr1, "FOS_chr1_mean", "MTS_chr1_mean", "chr1"), "FOS_chr1_mean", "Mean_TPM"),
  prepare_data(load_and_filter_data(gonad_files$chrx, "FOS_chrx_mean", "MTS_chrx_mean", "chrX"), "FOS_chrx_mean", "Mean_TPM")
)

# Generate color palette
colourCount <- 14
palette1 <- rev(colorRampPalette(brewer.pal(12, "Paired"))(colourCount))

# Function to generate boxplots
plot_boxplot <- function(data, value_col, title, palette_used = NULL) {
  data <- data %>% mutate(log2_value = log2(.data[[value_col]]))  # Avoid log(0)
  plot <- ggbetweenstats(
    data = data, x = gene_id, y = log2_value,
    xlab = "", ylab = "log2(TPM)",
    type = "np", plot.type = "box",
    p.adjust.method = "BH", conf.level = 0.99,
    pairwise.display = "significant", title = title,
    package = "ggsci", palette = "nrc_npg",
    boxplot.args = list(width = 0.5, alpha = 0.2, na.rm = TRUE),
    violin.args = list(width = 0.7, alpha = 0.2, na.rm = TRUE),
    centrality.label.args = list(size = 2, nudge_x = 0.4, segment.linetype = 4,
                                 min.segment.length = 0)
  )
  if (!is.null(palette_used)) {
    plot <- plot + scale_color_manual(name = "Chromosome", values = palette_used)
  }
  plot + theme(legend.position = 'none', legend.text = element_text(size = 10))
}

# Generate plots
p7 <- plot_boxplot(Male_combined, "Mean_TPM", "Male Gonads", palette1)
p8 <- plot_boxplot(Female_combined, "Mean_TPM", "Female Gonads", palette1)

p7$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)
p8$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)

print(p7)
print(p8)

# Compute Male/Female ratio
MF_ratio <- Male_combined %>%
  mutate(Ratio = Mean_TPM / Female_combined$Mean_TPM) %>%
  mutate(log2_Ratio = log2(Ratio))  # Avoid log(0)

# Plot Male/Female ratio
p9 <- ggbetweenstats(
  data = MF_ratio, x = gene_id, y = log2_Ratio,
  xlab = "", ylab = "log2(M/F TPM)",
  type = "np", plot.type = "box",
  p.adjust.method = "BH", conf.level = 0.99,
  pairwise.display = "significant", title = "Male to Female Ratio",
  package = "ggsci", palette = "nrc_npg",
  boxplot.args = list(width = 0.5, alpha = 0.2, na.rm = TRUE),
  violin.args = list(width = 0.7, alpha = 0.2, na.rm = TRUE),
  centrality.label.args = list(size = 2, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)
) + scale_color_manual(name = "Chromosome", values = palette1) + theme(legend.position = 'none')

p9$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)
print(p9)

# Combine plots using patchwork
combined_plot <- p1 + p2 + p5 + p3 + p4 + p6 + p7 + p8 + p9 +
  plot_annotation(tag_levels = list(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'))) &
  theme(plot.tag = element_text(face = 'bold'))
print(combined_plot)
