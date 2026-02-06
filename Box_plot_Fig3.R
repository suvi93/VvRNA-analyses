# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(ggstatsplot)
library(RColorBrewer)
library(patchwork)

setwd("/Users/octaviopalacios/Desktop/DC_paper/Figures/Figure_2")

# Function to load and filter data
load_and_filter_data <- function(file, female_col, male_col, chr_label) {
  df <- read.delim(file, header = TRUE, check.names = TRUE)
  print(paste("Loading file:", file))
  print("Column names detected:")
  print(colnames(df))
  
  if (!(female_col %in% colnames(df)) | !(male_col %in% colnames(df))) {
    stop(paste("Error: Column", female_col, "or", male_col, "not found in", file))
  }
  
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

Male_combined <- list()
Female_combined <- list()

body_parts <- list(
  "head" = list(
    "autosomes" = list("file" = "MF_head_ATPMmean.tsv", "female_col" = "FHS_A", "male_col" = "MHS_A"),
    "par" = list("file" = "P24XYMF_head_PARXrYTPMmean.tsv.new", "female_col" = "FHS_XrY", "male_col" = "MHS_XrY"),
    "slr" = list("file" = "P24XYMF_head_SLRXrYTPMmean.tsv.new", "female_col" = "FHS_XrY", "male_col" = "MHS_XrY"),
    "chrx" = list("file" = "MF_head_XlTPMmean.tsv", "female_col" = "FHS_Xl", "male_col" = "MHS_Xl")
  ),
  "legs" = list(
    "autosomes" = list("file" = "MF_legs_ATPMmean.tsv", "female_col" = "FLS_A", "male_col" = "MLS_A"),
    "par" = list("file" = "P24XYMF_legs_PARXrYTPMmean.tsv.new", "female_col" = "FLS_XrY", "male_col" = "MLS_XrY"),
    "slr" = list("file" = "P24XYMF_legs_SLRXrYTPMmean.tsv.new", "female_col" = "FLS_XrY", "male_col" = "MLS_XrY"),
    "chrx" = list("file" = "MF_legs_XlTPMmean.tsv", "female_col" = "FLS_Xl", "male_col" = "MLS_Xl")
  )
)

for (part in names(body_parts)) {
  datasets <- body_parts[[part]]
  
  Male_combined[[part]] <- bind_rows(
    prepare_data(load_and_filter_data(datasets$autosomes$file, datasets$autosomes$female_col, datasets$autosomes$male_col, "Autosomes"), datasets$autosomes$male_col, "Mean_TPM"),
    prepare_data(load_and_filter_data(datasets$par$file, datasets$par$female_col, datasets$par$male_col, "chrXR-Y PAR"), datasets$par$male_col, "Mean_TPM"),
    prepare_data(load_and_filter_data(datasets$slr$file, datasets$slr$female_col, datasets$slr$male_col, "chrXR-Y SLR"), datasets$slr$male_col, "Mean_TPM"),
    prepare_data(load_and_filter_data(datasets$chrx$file, datasets$chrx$female_col, datasets$chrx$male_col, "chrXL"), datasets$chrx$male_col, "Mean_TPM")
  )
  
  Female_combined[[part]] <- bind_rows(
    prepare_data(load_and_filter_data(datasets$autosomes$file, datasets$autosomes$female_col, datasets$autosomes$male_col, "Autosomes"), datasets$autosomes$female_col, "Mean_TPM"),
    prepare_data(load_and_filter_data(datasets$par$file, datasets$par$female_col, datasets$par$male_col, "chrXR-Y PAR"), datasets$par$female_col, "Mean_TPM"),
    prepare_data(load_and_filter_data(datasets$slr$file, datasets$slr$female_col, datasets$slr$male_col, "chrXR-Y SLR"), datasets$slr$female_col, "Mean_TPM"),
    prepare_data(load_and_filter_data(datasets$chrx$file, datasets$chrx$female_col, datasets$chrx$male_col, "chrXL"), datasets$chrx$female_col, "Mean_TPM")
  )
}

plot_boxplot <- function(data, value_col, title, palette_used = NULL) {
  data$gene_id <- factor(data$gene_id,
                         levels = c("Autosomes", "chrXR-Y PAR", "chrXR-Y SLR", "chrXL"))
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
  if (!is.null(palette_used)) {
    plot <- plot + scale_color_manual(name = "Chromosome", values = palette_used)
  }
  plot + theme(legend.position = 'none', legend.text = element_text(size = 10))
}

colourCount <- 14
palette1 <- rev(colorRampPalette(brewer.pal(12, "Paired"))(colourCount))

p1 <- plot_boxplot(Male_combined$head, "Mean_TPM", "Male Head", palette1)
p2 <- plot_boxplot(Female_combined$head, "Mean_TPM", "Female Head", palette1)
p3 <- plot_boxplot(Male_combined$legs, "Mean_TPM", "Male Leg", palette1)
p4 <- plot_boxplot(Female_combined$legs, "Mean_TPM", "Female Leg", palette1)

p1$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)
p2$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)
p3$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)
p4$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)

print(p1)
print(p2)
print(p3)
print(p4)

MF_ratio_head <- Male_combined$head %>% mutate(Ratio = Mean_TPM / Female_combined$head$Mean_TPM) %>% mutate(log2_Ratio = log2(Ratio))
MF_ratio_legs <- Male_combined$legs %>% mutate(Ratio = Mean_TPM / Female_combined$legs$Mean_TPM) %>% mutate(log2_Ratio = log2(Ratio))

MF_ratio_head$gene_id <- factor(MF_ratio_head$gene_id,
                                levels = c("Autosomes", "chrXR-Y PAR", "chrXR-Y SLR", "chrXL"))
MF_ratio_legs$gene_id <- factor(MF_ratio_legs$gene_id,
                                levels = c("Autosomes", "chrXR-Y PAR", "chrXR-Y SLR", "chrXL"))

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

gonad_files <- list(
  "autosomes" = "MF_gonads_ATPMmean.tsv",
  "par" = "P24XYMF_gonads_PARXrYTPMmean.tsv.new",
  "slr" = "P24XYMF_gonads_SLRXrYTPMmean.tsv.new",
  "chrx" = "MF_gonads_XlTPMmean.tsv"
)

Male_combined_gonads <- bind_rows(
  prepare_data(load_and_filter_data(gonad_files$autosomes, "FOS_A", "MTS_A", "Autosomes"), "MTS_A", "Mean_TPM"),
  prepare_data(load_and_filter_data(gonad_files$par, "FOS_XrY", "MTS_XrY", "chrXR-Y PAR"), "MTS_XrY", "Mean_TPM"),
  prepare_data(load_and_filter_data(gonad_files$slr, "FOS_XrY", "MTS_XrY", "chrXR-Y SLR"), "MTS_XrY", "Mean_TPM"),
  prepare_data(load_and_filter_data(gonad_files$chrx, "FOS_Xl", "MTS_Xl", "chrXL"), "MTS_Xl", "Mean_TPM")
)

Female_combined_gonads <- bind_rows(
  prepare_data(load_and_filter_data(gonad_files$autosomes, "FOS_A", "MTS_A", "Autosomes"), "FOS_A", "Mean_TPM"),
  prepare_data(load_and_filter_data(gonad_files$par, "FOS_XrY", "MTS_XrY", "chrXR-Y PAR"), "FOS_XrY", "Mean_TPM"),
  prepare_data(load_and_filter_data(gonad_files$slr, "FOS_XrY", "MTS_XrY", "chrXR-Y SLR"), "FOS_XrY", "Mean_TPM"),
  prepare_data(load_and_filter_data(gonad_files$chrx, "FOS_Xl", "MTS_Xl", "chrXL"), "FOS_Xl", "Mean_TPM")
)

p7 <- plot_boxplot(Male_combined_gonads, "Mean_TPM", "Male Gonad", palette1)
p8 <- plot_boxplot(Female_combined_gonads, "Mean_TPM", "Female Gonad", palette1)

p7$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)
p8$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)

print(p7)
print(p8)

MF_ratio_gonads <- Male_combined_gonads %>%
  mutate(Ratio = Mean_TPM / Female_combined_gonads$Mean_TPM) %>%
  mutate(log2_Ratio = log2(Ratio))

MF_ratio_gonads$gene_id <- factor(MF_ratio_gonads$gene_id,
                                  levels = c("Autosomes", "chrXR-Y PAR", "chrXR-Y SLR", "chrXL"))

p9 <- ggbetweenstats(
  data = MF_ratio_gonads, x = gene_id, y = log2_Ratio,
  xlab = "", ylab = "log2(M/F TPM)",
  type = "np", plot.type = "box",
  p.adjust.method = "BH", conf.level = 0.99,
  pairwise.display = "significant", title = "Male to Female Ratio - Gonad",
  package = "ggsci", palette = "nrc_npg",
  boxplot.args = list(width = 0.5, alpha = 0.2, na.rm = TRUE),
  violin.args = list(width = 0.7, alpha = 0.2, na.rm = TRUE),
  centrality.label.args = list(size = 2, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)
  
) + scale_color_manual(name = "Chromosome", values = palette1) + theme(legend.position = 'none')

# Modify the jitter layer inside the ggplot object
# Identify which layer uses position_jitter and change its width
p9$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)

print(p9)

combined_plot <- p1 + p2 + p5 + p3 + p4 + p6 + p7 + p8 + p9 +
  plot_annotation(tag_levels = list(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'))) &
  theme(plot.tag = element_text(face = 'bold'))
print(combined_plot)