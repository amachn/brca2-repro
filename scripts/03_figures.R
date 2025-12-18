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
# - 1.0 ~ Setup (Data, Style, Funcs)
df_tcga <- as.data.frame(colData(se_tcga)) |>
  mutate(
    BRCA2_group = factor(BRCA2_group, levels = c("Low", "High")),
    ajcc_stage_simple = case_when(
      str_detect(ajcc_stage, "^Stage IV") ~ "IV",
      str_detect(ajcc_stage, "^Stage III") ~ "III",
      str_detect(ajcc_stage, "^Stage II") ~ "II",
      str_detect(ajcc_stage, "^Stage I") ~ "I",
      TRUE ~ NA_character_
    ),
    ajcc_stage_simple = factor(
      ajcc_stage_simple, levels = c("I", "II", "III", "IV")
    )
  )
df_gse <- as.data.frame(colData(se_gse)) |>
  mutate(BRCA2_group = factor(BRCA2_group, levels = c("Low", "High")))

p_box <- function(df, y, x="BRCA2_group", test="wilcox.test",
                  title=NULL, ylab=NULL, show_title=TRUE, show_y=FALSE) {
  ggplot(df, aes(.data[[x]], .data[[y]])) +
    geom_boxplot(outlier.size = 0.6) +
    stat_compare_means(
      method = test,
      label = "p.format",
      label.x.npc = "center",
      label.y.npc = "bottom",
      hjust = 0.5
    ) +
    labs(title = if(show_title) title else NULL,
         y = if(show_y) ylab else NULL,
         x = NULL) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, vjust = 1, size = 10),
      axis.title.y = element_text(size = 9),
      plot.margin = margin(14, 10, 10, 10)
    )
}

cohort_tag <- function(lbl) {
  ggplot() +
    annotate("text", x = 0, y = 0.5, label = lbl, fontface = "bold", size = 4) +
    theme_void()
}

# - 1.1 ~ Figure 2
# - 1.1.1 ~ Top Row (TCGA): AJCC, Subtype, MKI67
p2_tcga_ajcc <- p_box(
  df_tcga |> filter(!is.na(ajcc_stage_simple)),
  y = "BRCA2", x = "ajcc_stage_simple",
  test = "kruskal.test", title = "AJCC Stage",
  ylab = "BRCA2 expression", show_y = TRUE
)

p2_tcga_subtype <- p_box(
  df_tcga |> filter(!is.na(Subtype)),
  y = "BRCA2", x = "Subtype",
  test = "kruskal.test", title = "Subtype"
)

p2_tcga_mki67 <- p_box(
  df_tcga, y = "MKI67",
  title = "MKI67", ylab = "MKI67 expression", show_y = TRUE
)

# - 1.1.2 ~ Bottom Row (GSE): NHG, Subtype, MKI67
p2_gse_nhg <- p_box(
  df_gse |> filter(nhg %in% c("G1", "G2", "G3")),
  y = "BRCA2", x = "nhg",
  test = "kruskal.test", title = "NHG",
  ylab = "BRCA2 expression", show_y = TRUE
)

p2_gse_subtype <- p_box(
  df_gse |> filter(!is.na(Subtype)),
  y = "BRCA2", x = "Subtype",
  test = "kruskal.test"
)

p2_gse_mki67 <- p_box(df_gse, y = "MKI67")

# - 1.1.3 ~ stacking
fig2 <- (cohort_tag("TCGA-BRCA") | p2_tcga_ajcc | p2_tcga_subtype | p2_tcga_mki67) /
        (cohort_tag("GSE96058")  | p2_gse_nhg   | p2_gse_subtype  | p2_gse_mki67)

# - 1.2 ~ Figure 3
# - 1.2.1 ~ 3a. E2F1/RB1/PALB2/PARP1 expr
p3a_tcga <- list(
  p_box(df_tcga, "E2F1", title = "E2F1",
        ylab = "mRNA expression", show_y = TRUE),
  p_box(df_tcga, "RB1", title = "RB1"),
  p_box(df_tcga, "PALB2", title = "PALB2"),
  p_box(df_tcga, "PARP1", title = "PARP1")
)

p3a_gse <- list(
  p_box(df_gse, "E2F1", show_title = FALSE,
        ylab = "mRNA expression", show_y = TRUE),
  p_box(df_gse, "RB1", show_title = FALSE),
  p_box(df_gse, "PALB2", show_title = FALSE),
  p_box(df_gse, "PARP1", show_title = FALSE)
)

# - 1.2.2 ~ 3b. E2F/G2M/P53
p3b_tcga <- list(
  p_box(df_tcga, "E2F", title = "E2F targets",
        ylab = "Score", show_y = TRUE),
  p_box(df_tcga, "G2M", title = "G2M checkpoint"),
  p_box(df_tcga, "P53", title = "P53 pathway")
)

p3b_gse <- list(
  p_box(df_gse, "E2F", show_title = FALSE,
        ylab = "Score", show_y = TRUE),
  p_box(df_gse, "G2M", show_title = FALSE),
  p_box(df_gse, "P53", show_title = FALSE)
)

# - 1.2.3 ~ 3c. ESR1
p3c_tcga <- p_box(df_tcga, "ESR1", title = "ESR1")
p3c_gse  <- p_box(df_gse, "ESR1", show_title = FALSE)

# - 1.2.4 ~ 3d. ER Response
p3d_tcga <- list(
  p_box(df_tcga, "ER_early", title = "ER response early"),
  p_box(df_tcga, "ER_late", title = "ER response late")
)

p3d_gse <- list(
  p_box(df_gse, "ER_early", show_title = FALSE),
  p_box(df_gse, "ER_late", show_title = FALSE)
)

# - 1.2.5 ~ stacking
row1 <- cohort_tag("TCGA-BRCA") | p3a_tcga[[1]] | p3a_tcga[[2]] | p3a_tcga[[3]] | p3a_tcga[[4]] | p3c_tcga
row2 <- cohort_tag("GSE96058")  | p3a_gse[[1]]  | p3a_gse[[2]]  | p3a_gse[[3]]  | p3a_gse[[4]]  | p3c_gse
row3 <- cohort_tag("TCGA-BRCA") | p3b_tcga[[1]] | p3b_tcga[[2]] | p3b_tcga[[3]] | p3d_tcga[[1]] | p3d_tcga[[2]]
row4 <- cohort_tag("GSE96058")  | p3b_gse[[1]]  | p3b_gse[[2]]  | p3b_gse[[3]]  | p3d_gse[[1]]  | p3d_gse[[2]]
fig3 <- (row1 / row2 / row3 / row4)
