library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# P24X0 data
data1 <- data.frame(
  Tissue = c("Head", "Legs", "Gonads"),
  MBG_X = c(0, 0, 14),
  MBG_1 = c(4, 3, 99),
  MBG_A = c(2, 2, 767),
  UBG_X = c(1226, 1126, 1014),
  UBG_1 = c(1158, 1062, 977),
  UBG_A = c(10132, 9325, 8560),
  FBG_X = c(2, 3, 144),
  FBG_1 = c(0, 0, 55),
  FBG_A = c(1, 5, 601)
)

data_long1 <- data1 %>%
  pivot_longer(cols = -Tissue, names_to = c("Type", "Group"), names_sep = "_") %>%
  group_by(Tissue, Group) %>% # Group by Tissue and Group to calculate proportions correctly
  mutate(Proportion = value / sum(value)) %>%
  ungroup() %>%
  mutate(Type = factor(Type, levels = c("MBG", "UBG", "FBG")),
         Group = factor(Group, levels = c("A", "1", "X"))) %>%
  arrange(Tissue, Group)

p1 <- ggplot(data_long1, aes(x = Group, y = Proportion, fill = Type)) +
        geom_bar(stat = "identity", position = "stack", color = "black") +
        scale_fill_manual(values = c("MBG" = "orange", "UBG" = "gray", "FBG" = "darkgreen"),
                          labels = c("MBG" = "Male-biased genes", "UBG" = "Unbiased genes", "FBG" = "Female-biased genes")) +
        facet_wrap(~ Tissue, scales = "free_x", ncol = 3) +
        theme_minimal() +
        theme(
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line.x = element_line(color = "black"),
              axis.line.y = element_line(color = "black"),
              axis.ticks = element_line(color = "black"),
              legend.position = "none",
              strip.text = element_text(size = 12, face = "bold"),
              axis.text.x = element_text(size = 8, face = "bold", color = "black"),
              axis.text.y = element_text(size = 8, color = "black"),
              plot.title = element_text(hjust = 0.5, face = "bold")
             ) +
        labs(x = "", y = "Proportion", fill = "", title = "Proportion of sex-biased genes in P24X0" ) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
        scale_x_discrete(labels=c("Autosomes", "chr1", "chrX"))

# P24XY data
data2 <- data.frame(
  Tissue = c("Head", "Legs", "Gonads"),
  MBG_XL = c(0, 0, 16),
  MBG_XRYPAR = c(0, 2, 92),
  MBG_XRYSLR = c(0, 0, 18),
  MBG_A = c(0, 2, 893),
  UBG_XL = c(1045, 965, 815),
  UBG_XRYPAR = c(926, 816, 769),
  UBG_XRYSLR = c(179, 148, 148),
  UBG_A = c(9719, 8628, 8211),
  FBG_XL = c(2, 2, 198),
  FBG_XRYPAR = c(0, 0, 59),
  FBG_XRYSLR = c(0, 0, 11),
  FBG_A = c(6, 9, 723)
)

data_long2 <- data2 %>%
  pivot_longer(
    cols = -Tissue,
    names_to = c("Type", "Group"),
    names_pattern = "([A-Z]+)_([A-Z_]+)"
  ) %>%
  group_by(Tissue, Group) %>%
  mutate(Proportion = value / sum(value)) %>%
  ungroup() %>%
  mutate(
    Type = factor(Type, levels = c("MBG", "UBG", "FBG")),
    Group = factor(Group, levels = c("A", "XRYPAR", "XRYSLR", "XL"))
  ) %>%
  arrange(Tissue, Group)

p2 <- ggplot(data_long2, aes(x = Group, y = Proportion, fill = Type)) +
        geom_bar(stat = "identity", position = "stack", color = "black") +
        scale_fill_manual(values = c("MBG" = "orange", "UBG" = "gray", "FBG" = "darkgreen"),
        facet_wrap(~ Tissue, scales = "free_x", ncol = 4) +
        theme_minimal() +
        theme(
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line.x = element_line(color = "black"),
              axis.line.y = element_line(color = "black"),
              axis.ticks = element_line(color = "black"),
              legend.position = "bottom",
              strip.text = element_text(size = 12, face = "bold"),
              axis.text.x = element_text(size = 5, face = "bold", color = "black"),
              axis.text.y = element_text(size = 8, color = "black"),
              plot.title = element_text(hjust = 0.5, face = "bold")
             ) +
        labs(x = "", y = "Proportion", fill = "", title = "Proportion of sex-biased genes in P24XY" ) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
        scale_x_discrete(labels=c("Autosomes", "chrXR-Y PAR", "chrXR-Y SLR", "chrXL"))

pdf("sbg_xavstackedbar.pdf")
grid.arrange(p1, p2, nrow = 2)
dev.off()
