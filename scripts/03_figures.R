# ~~~ 03_figures.R ~~~
# - 0 ~ Libraries + Data
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ggpubr)
  library(tidyverse)
})

load("data/processed/SE_tcga_gse_anno.RData")

# - 1 ~ FIGURES
# - 1.0 ~ Setup + Helper Functions
df_tcga <- as.data.frame(colData(se_tcga))
df_gse  <- as.data.frame(colData(se_gse))
df_tcga$BRCA2_group <- factor(df_tcga$BRCA2_group, levels = c("Low", "High"))
df_gse$BRCA2_group  <- factor(df_gse$BRCA2_group,  levels = c("Low", "High"))

p_box <- function(df, x, y, title = NULL, test = c("wilcox.test", "kruskal.test")) {
  test <- match.arg(test)
  ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_boxplot() +
    stat_compare_means(method = test) +
    theme_bw() +
    labs(title = title, x = NULL, y = NULL)
}

# - 1.1 ~ Figure 1

# - 1.2 ~ Figure 2
# - 1.2.1 ~ 2a.

# - 1.2.2 ~ 2b.
p2b1 <- p_box(df_tcga, "BRCA2_group", "MKI67",
              "TCGA: MKI67 by BRCA2 group", test = "wilcox.test")
p2b2 <- p_box(df_gse, "BRCA2_group", "MKI67",
              "GSE: MKI67 by BRCA2 group", test = "wilcox.test")

# - 1.3 ~ Figure 3
