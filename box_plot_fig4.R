# ------------------------------------------------------------------
# box_plot_fig4.R
#
# Generates Figure 4: combined boxplot panels comparing expression
# patterns across tissues and chromosomal regions between P24X0 and
# P24XY races.
#
# Input:  per-tissue/region TPM mean tables (derived from deseq2.R /
#         featureCounts output)
# Output: combined 3-panel figure (Figure 4 of the paper)
#
# Part of: VvRNA-analyses pipeline (step 3 — visualization)
# Author:  Octavio Palacios, Suvratha Jayaprasad
# Paper:   Jayaprasad et al. 2026, Genome Biology and Evolution,
#          https://doi.org/10.1093/gbe/evag026
# License: MIT — see LICENSE
# ------------------------------------------------------------------

library(readr)
library(dplyr)
library(ggplot2)
library(ggstatsplot)
library(RColorBrewer)
library(patchwork)

# ====================================================================
# GONADS
# ====================================================================

# Load existing autosomes data (P24X0)
P24XOMF <- read_delim("P24XOMF_gonads_autosomes_TPMmean.tsv.new",
                      delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_autosomes <- P24XOMF %>%
  filter(FOS_autosomes_mean > 1 & MTS_autosomes_mean > 1) %>%
  mutate(log2_ratio = log2(MTS_autosomes_mean / FOS_autosomes_mean)) %>%
  select(gene_id, log2_ratio) %>%
  mutate(group = "P24X0 autosomes")

# Load XY autosomes data (P24XY)
P24XYMF_gonads_ATPMmean_tsv <- read_delim("P24XYMF_gonads_ATPMmean.tsv.new",
                                          delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_XY <- P24XYMF_gonads_ATPMmean_tsv %>%
  filter(FOS_A > 1 & MTS_A > 1) %>%
  mutate(log2_ratio = log2(MTS_A / FOS_A)) %>%
  select(gene_ids, log2_ratio) %>%
  rename(gene_id = gene_ids) %>%
  mutate(group = "P24XY autosomes")

# Load chr1 data (P24X0)
P24XOMF_gonads_chr1 <- read_delim("P24XOMF_gonads_chr1_TPMmean.tsv.new",
                                  delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_chr1 <- P24XOMF_gonads_chr1 %>%
  filter(FOS_chr1_mean > 1 & MTS_chr1_mean > 1) %>%
  mutate(log2_ratio = log2(MTS_chr1_mean / FOS_chr1_mean)) %>%
  select(gene_id, log2_ratio) %>%
  mutate(group = "P24X0 chr1")

# Load XrY data (P24XY)
P24XYMF_gonads_XrYPAR <- read_delim("P24XYMF_gonads_XrYPARTPMmean.tsv.new",
                                 delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_XrYPAR <- P24XYMF_gonads_XrYPAR %>%
  filter(FOS_XrY > 1 & MTS_XrY > 1) %>%
  mutate(log2_ratio = log2(MTS_XrY / FOS_XrY)) %>%
  select(gene_ids, log2_ratio) %>%
  rename(gene_id = gene_ids) %>%
  mutate(group = "P24XY chr1 (XR-Y PAR)")

P24XYMF_gonads_XrYSLR <- read_delim("P24XYMF_gonads_XrYSLRTPMmean.tsv.new",
                                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_XrYSLR <- P24XYMF_gonads_XrYSLR %>%
  filter(FOS_XrY > 1 & MTS_XrY > 1) %>%
  mutate(log2_ratio = log2(MTS_XrY / FOS_XrY)) %>%
  select(gene_ids, log2_ratio) %>%
  rename(gene_id = gene_ids) %>%
  mutate(group = "P24XY chr1 (XR-Y SLR)")

# Load P24X0 chrX data
P24XOMF_gonads_chrx <- read_delim("P24XOMF_gonads_chrx_TPMmean.tsv.new",
                                  delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_chrx <- P24XOMF_gonads_chrx %>%
  rename(gene_id = gene_id) %>%
  filter(FOS_chrx_mean > 1 & MTS_chrx_mean > 1) %>%
  mutate(log2_ratio = log2(MTS_chrx_mean / FOS_chrx_mean)) %>%
  select(gene_id, log2_ratio) %>%
  mutate(group = "P24X0 chrX")

# Load P24XY chrXl data
P24XYMF_gonads_Xl <- read_delim("P24XYMF_gonads_XlTPMmean.tsv.new",
                                delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_Xl <- P24XYMF_gonads_Xl %>%
  filter(FOS_Xl > 1 & MTS_Xl > 1) %>%
  mutate(log2_ratio = log2(MTS_Xl / FOS_Xl)) %>%  # Ensure FOS_Xl exists
  select(gene_ids, log2_ratio) %>%
  rename(gene_id = gene_ids) %>%
  mutate(group = "P24XY chrX (XL)")

# Combine all datasets
combined_data <- bind_rows(
  filtered_autosomes,
  filtered_XY,
  filtered_chr1,
  filtered_XrYPAR,
  filtered_XrYSLR,
  filtered_chrx,
  filtered_Xl
) %>%
  mutate(group = factor(group, levels = c(
    "P24X0 autosomes",
    "P24XY autosomes",
    "P24X0 chr1",
    "P24XY chr1 (XR-Y PAR)",
    "P24XY chr1 (XR-Y SLR)",
    "P24X0 chrX",
    "P24XY chrX (XL)"
  )))

# Plot combined data with specified order
plot_combined_1 <- ggbetweenstats(
  data = combined_data,
  x = group,
  y = log2_ratio,
  xlab = "",
  ylab = "log2(Male/Female TPM)",
  type = "np",
  plot.type = "box",
  p.adjust.method = "BH",
  conf.level = 0.99,
  pairwise.display = "s",
  title = "Male to Female Ratio - Gonad",
  package = "ggsci",
  palette = "nrc_npg",
  boxplot.args = list(width = 0.5, alpha = 0.2, na.rm = TRUE),
  violin.args = list(width = 0.7, alpha = 0.2, na.rm = TRUE),
  centrality.label.args = list(size = 2, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)
)

# Modify the jitter layer inside the ggplot object
# Identify which layer uses position_jitter and change its width
plot_combined_1$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)

plot_combined_1

# ====================================================================
# HEAD
# ====================================================================

# Load existing head autosomes data (P24X0)
P24XOMF_head_autosomes_TPMmean.tsv <- read_delim("P24XOMF_head_autosomes_TPMmean.tsv.new",
                                                 delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_head_autosomes <- P24XOMF_head_autosomes_TPMmean.tsv %>%
  filter(FHS_autosomes_mean > 1 & MHS_autosomes_mean > 1) %>%
  mutate(log2_ratio = log2(MHS_autosomes_mean / FHS_autosomes_mean)) %>%
  select(gene_id, log2_ratio) %>%
  mutate(group = "P24X0 autosomes")

# Load XY head autosomes data (P24XY)
P24XYMF_head_ATPMmean.tsv <- read_delim("P24XYMF_head_ATPMmean.tsv.new",
                                        delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_head_XY <- P24XYMF_head_ATPMmean.tsv %>%
  filter(FHS_A > 1 & MHS_A > 1) %>%
  mutate(log2_ratio = log2(MHS_A / FHS_A)) %>%
  select(gene_ids, log2_ratio) %>%
  rename(gene_id = gene_ids) %>%
  mutate(group = "P24XY autosomes")

# Load chr1 head data (P24X0)
P24XOMMF_head_chr1_TPMmean.tsv <- read_delim("P24XOMMF_head_chr1_TPMmean.tsv.new",
                                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_head_chr1 <- P24XOMMF_head_chr1_TPMmean.tsv %>%
  filter(FHS_chr1_mean > 1 & MHS_chr1_mean > 1) %>%
  mutate(log2_ratio = log2(MHS_chr1_mean / FHS_chr1_mean)) %>%
  select(gene_id, log2_ratio) %>%
  mutate(group = "P24X0 chr1")

# Load XrY head data (P24XY)
P24XYMF_head_XrYPARTPMmean.tsv <- read_delim("P24XYMF_head_XrYPARTPMmean.tsv.new",
                                          delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_head_XrYPAR <- P24XYMF_head_XrYPARTPMmean.tsv %>%
  filter(FHS_XrY > 1 & MHS_XrY > 1) %>%
  mutate(log2_ratio = log2(MHS_XrY / FHS_XrY)) %>%
  select(gene_ids, log2_ratio) %>%
  rename(gene_id = gene_ids) %>%
  mutate(group = "P24XY chr1 (XR-Y PAR)")

P24XYMF_head_XrYSLRTPMmean.tsv <- read_delim("P24XYMF_head_XrYSLRTPMmean.tsv.new",
                                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_head_XrYSLR <- P24XYMF_head_XrYSLRTPMmean.tsv %>%
  filter(FHS_XrY > 1 & MHS_XrY > 1) %>%
  mutate(log2_ratio = log2(MHS_XrY / FHS_XrY)) %>%
  select(gene_ids, log2_ratio) %>%
  rename(gene_id = gene_ids) %>%
  mutate(group = "P24XY chr1 (XR-Y SLR)")

# Load P24X0 chrX head data
P24XOMF_head_chrx_TPMmean.tsv <- read_delim("P24XOMF_head_chrx_TPMmean.tsv.new",
                                            delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_head_chrx <- P24XOMF_head_chrx_TPMmean.tsv %>%
  rename(gene_id = gene_id) %>%
  filter(FHS_chrx_mean > 1 & MHS_chrx_mean > 1) %>%
  mutate(log2_ratio = log2(MHS_chrx_mean / FHS_chrx_mean)) %>%
  select(gene_id, log2_ratio) %>%
  mutate(group = "P24X0 chrX")

# Load P24XY chrXl head data
P24XYMF_head_XlTPMmean.tsv <- read_delim("P24XYMF_head_XlTPMmean.tsv.new",
                                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_head_Xl <- P24XYMF_head_XlTPMmean.tsv %>%
  filter(FHS_Xl > 1 & MHS_Xl > 1) %>%
  mutate(log2_ratio = log2(MHS_Xl / FHS_Xl)) %>%  # Ensure FOS_Xl exists
  select(gene_ids, log2_ratio) %>%
  rename(gene_id = gene_ids) %>%
  mutate(group = "P24XY chrX (XL)")

# Combine all datasets
combined_data_2 <- bind_rows(
  filtered_head_autosomes,
  filtered_head_XY,
  filtered_head_chr1,
  filtered_head_XrYPAR,
  filtered_head_XrYSLR,
  filtered_head_chrx,
  filtered_head_Xl
) %>%
  mutate(group = factor(group, levels = c(
    "P24X0 autosomes",
    "P24XY autosomes",
    "P24X0 chr1",
    "P24XY chr1 (XR-Y PAR)",
    "P24XY chr1 (XR-Y SLR)",
    "P24X0 chrX",
    "P24XY chrX (XL)"
  )))

# Plot combined data with specified order
plot_combined_2 <- ggbetweenstats(
  data = combined_data_2,
  x = group,
  y = log2_ratio,
  xlab = "",
  ylab = "log2(Male/Female TPM)",
  type = "np",
  plot.type = "box",
  p.adjust.method = "BH",
  conf.level = 0.99,
  pairwise.display = "s",
  title = "Male to Female Ratio - Head",
  package = "ggsci",
  palette = "nrc_npg",
  boxplot.args = list(width = 0.5, alpha = 0.2, na.rm = TRUE),
  violin.args = list(width = 0.7, alpha = 0.2, na.rm = TRUE),
  centrality.label.args = list(size = 2, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)
)

# Modify the jitter layer inside the ggplot object
# Identify which layer uses position_jitter and change its width
plot_combined_2$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)

plot_combined_2

# ====================================================================
# LEGS
# ====================================================================

# Load existing leg autosomes data (P24X0)
P24XOMF_legs_autosomes_TPMmean.tsv <- read_delim("P24XOMF_legs_autosomes_TPMmean.tsv.new", 
                                                 delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_leg_autosomes <- P24XOMF_legs_autosomes_TPMmean.tsv %>%
  filter(FLS_autosomes_mean > 1 & MLS_autosomes_mean > 1) %>%
  mutate(log2_ratio = log2(MLS_autosomes_mean / FLS_autosomes_mean)) %>%
  select(gene_id, log2_ratio) %>%
  mutate(group = "P24X0 autosomes")

# Load XY leg autosomes data (P24XY)
P24XYMF_legs_ATPMmean.tsv <- read_delim("P24XYMF_legs_ATPMmean.tsv.new",
                                        delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_leg_XY <- P24XYMF_legs_ATPMmean.tsv %>%
  filter(FLS_A > 1 & MLS_A > 1) %>%
  mutate(log2_ratio = log2(MLS_A / FLS_A)) %>%
  select(gene_ids, log2_ratio) %>%
  rename(gene_id = gene_ids) %>%
  mutate(group = "P24XY autosomes")

# Load chr1 leg data (P24X0)
P24XOMF_legs_chr1_TPMmean.tsv <- read_delim("P24XOMF_legs_chr1_TPMmean.tsv.new",
                                            delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_leg_chr1 <- P24XOMF_legs_chr1_TPMmean.tsv %>%
  filter(FLS_chr1_mean > 1 & MLS_chr1_mean > 1) %>%
  mutate(log2_ratio = log2(MLS_chr1_mean / FLS_chr1_mean)) %>%
  select(gene_id, log2_ratio) %>%
  mutate(group = "P24X0 chr1")

# Load XrY head data (P24XY)
P24XYMF_legs_XrYPARTPMmean.tsv <- read_delim("P24XYMF_legs_XrYPARTPMmean.tsv.new",
                                          delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_leg_XrYPAR <- P24XYMF_legs_XrYPARTPMmean.tsv %>%
  filter(FLS_XrY > 1 & MLS_XrY > 1) %>%
  mutate(log2_ratio = log2(MLS_XrY / FLS_XrY)) %>%
  select(gene_ids, log2_ratio) %>%
  rename(gene_id = gene_ids) %>%
  mutate(group = "P24XY chr1 (XR-Y PAR)")

P24XYMF_legs_XrYSLRTPMmean.tsv <- read_delim("P24XYMF_legs_XrYSLRTPMmean.tsv.new",
                                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_leg_XrYSLR <- P24XYMF_legs_XrYSLRTPMmean.tsv %>%
  filter(FLS_XrY > 1 & MLS_XrY > 1) %>%
  mutate(log2_ratio = log2(MLS_XrY / FLS_XrY)) %>%
  select(gene_ids, log2_ratio) %>%
  rename(gene_id = gene_ids) %>%
  mutate(group = "P24XY chr1 (XR-Y SLR)")

# Load P24X0 chrX leg data
P24XOMF_legs_chrx_TPMmean.tsv <- read_delim("P24XOMF_legs_chrx_TPMmean.tsv.new",
                                            delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_leg_chrx <- P24XOMF_legs_chrx_TPMmean.tsv %>%
  rename(gene_id = gene_id) %>%
  filter(FLS_chrx_mean > 1 & MLS_chrx_mean > 1) %>%
  mutate(log2_ratio = log2(MLS_chrx_mean / FLS_chrx_mean)) %>%
  select(gene_id, log2_ratio) %>%
  mutate(group = "P24X0 chrX")

# Load P24XY chrXl leg data
P24XYMF_legs_XlTPMmean.tsv <- read_delim("P24XYMF_legs_XlTPMmean.tsv.new",
                                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)

filtered_leg_Xl <- P24XYMF_legs_XlTPMmean.tsv %>%
  filter(FLS_Xl > 1 & MLS_Xl > 1) %>%
  mutate(log2_ratio = log2(MLS_Xl / FLS_Xl)) %>%  # Ensure FOS_Xl exists
  select(gene_ids, log2_ratio) %>%
  rename(gene_id = gene_ids) %>%
  mutate(group = "P24XY chrX (XL)")

# Combine all datasets
combined_data_3 <- bind_rows(
  filtered_leg_autosomes,
  filtered_leg_XY,
  filtered_leg_chr1,
  filtered_leg_XrYPAR,
  filtered_leg_XrYSLR,
  filtered_leg_chrx,
  filtered_leg_Xl
) %>%
  mutate(group = factor(group, levels = c(
    "P24X0 autosomes",
    "P24XY autosomes",
    "P24X0 chr1",
    "P24XY chr1 (XR-Y PAR)",
    "P24XY chr1 (XR-Y SLR)",
    "P24X0 chrX",
    "P24XY chrX (XL)"
  )))

# Plot combined data with specified order
plot_combined_3 <- ggbetweenstats(
  data = combined_data_3,
  x = group,
  y = log2_ratio,
  xlab = "",
  ylab = "log2(Male/Female TPM)",
  type = "np",
  plot.type = "box",
  p.adjust.method = "BH",
  conf.level = 0.99,
  pairwise.display = "s",
  title = "Male to Female Ratio - Leg",
  package = "ggsci",
  palette = "nrc_npg",
  boxplot.args = list(width = 0.5, alpha = 0.2, na.rm = TRUE),
  violin.args = list(width = 0.7, alpha = 0.2, na.rm = TRUE),
  centrality.label.args = list(size = 2, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)
)

# Modify the jitter layer inside the ggplot object
# Identify which layer uses position_jitter and change its width
plot_combined_3$layers[[1]]$position <- position_jitter(width = 0.2, height = 0)

plot_combined_3

# ====================================================================
# COMBINE ALL PANELS INTO FIGURE 4
# ====================================================================

# Put the legend of x-axis 45-degrees
plot_combined_1 <- plot_combined_1 + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

plot_combined_2 <- plot_combined_2 + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

plot_combined_3 <- plot_combined_3 + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Combine all plots into one figure using patchwork
combined_final_plot <- plot_combined_2 + plot_combined_3 + plot_combined_1 +
  plot_annotation(tag_levels = list(c('a', 'b', 'c'))) &
  theme(plot.tag = element_text(face = 'bold'))
print(combined_final_plot)
