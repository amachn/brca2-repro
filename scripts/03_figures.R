# ~~~ 03_figures.R ~~~
# - 0 ~ Libraries + Data
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ggpubr)
  library(patchwork)
  library(tidyverse)
})

load("data/processed/SE_tcga_gse_anno.RData")

# - 1 ~ FIGURES
# - 1.0.a ~ Data Setup
df_tcga <- as.data.frame(colData(se_tcga)) |>
  mutate(BRCA2_group = factor(BRCA2_group, levels = c("Low", "High")))
df_gse  <- as.data.frame(colData(se_gse)) |>
  mutate(BRCA2_group = factor(BRCA2_group, levels = c("Low", "High")))

# - 1.0.b ~ Plot Style + Funcs
theme_paper <- function() {
  theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "italic", size = 10),
      axis.title = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
}

p_box <- function(df, y, title) {
  ggplot(df, aes(BRCA2_group, .data[[y]])) +
    geom_boxplot(outlier.size = 0.6) +
    stat_compare_means(method = "wilcox.test", label = "p.format", label.y.npc = "top", size = 3) +
    labs(title = title) +
    theme_paper()
}

cohort_tag <- function(lbl) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = lbl, angle = 90, fontface = "bold", size = 4) +
    theme_void() +
    theme(
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      plot.margin = margin(4, 4, 4, 4)
    )
}

# - 1.3 ~ Figure 3
# - 1.3.1 ~ 3a. E2F1/RB1/PALB2/PARP1 expr
p3a_tcga <- list(
  p_box(df_tcga, "E2F1", "E2F1"),
  p_box(df_tcga, "RB1", "RB1"),
  p_box(df_tcga, "PALB2", "PALB2"),
  p_box(df_tcga, "PARP1", "PARP1")
)

p3a_gse <- list(
  p_box(df_gse, "E2F1", "E2F1"),
  p_box(df_gse, "RB1", "RB1"),
  p_box(df_gse, "PALB2", "PALB2"),
  p_box(df_gse, "PARP1", "PARP1")
)

# - 1.3.2 ~ 3b. E2F/G2M/P53
p3b_tcga <- list(
  p_box(df_tcga, "E2F", "E2F targets"),
  p_box(df_tcga, "G2M", "G2M checkpoint"),
  p_box(df_tcga, "P53", "P53 pathway")
)

p3b_gse <- list(
  p_box(df_gse, "E2F", "E2F targets"),
  p_box(df_gse, "G2M", "G2M checkpoint"),
  p_box(df_gse, "P53", "P53 pathway")
)

# - 1.3.3 ~ 3c. ESR1
p3c_tcga <- p_box(df_tcga, "ESR1", "ESR1")
p3c_gse  <- p_box(df_gse, "ESR1", "ESR1")

# - 1.3.4 ~ 3d. ER
p3d_tcga <- list(
  p_box(df_tcga, "ER_early", "ER response early"),
  p_box(df_tcga, "ER_late", "ER response late")
)

p3d_gse <- list(
  p_box(df_gse, "ER_early", "ER response early"),
  p_box(df_gse, "ER_late", "ER response late")
)

# - 1.3.5 ~ stacking
row1 <- cohort_tag("TCGA")     | p3a_tcga[[1]] | p3a_tcga[[2]] | p3a_tcga[[3]] | p3a_tcga[[4]] | p3c_tcga
row2 <- cohort_tag("GSE96058") | p3a_gse[[1]]  | p3a_gse[[2]]  | p3a_gse[[3]]  | p3a_gse[[4]]  | p3c_gse
row3 <- cohort_tag("TCGA")     | p3b_tcga[[1]] | p3b_tcga[[2]] | p3b_tcga[[3]] | p3d_tcga[[1]] | p3d_tcga[[2]]
row4 <- cohort_tag("GSE96058") | p3b_gse[[1]]  | p3b_gse[[2]]  | p3b_gse[[3]]  | p3d_gse[[1]]  | p3d_gse[[2]]

fig3 <- (row1 / row2 / row3 / row4) +
  plot_layout(widths = c(0.45, 1, 1, 1, 1, 1))
fig3
