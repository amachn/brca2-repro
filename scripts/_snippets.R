# ~~~ _snippets.R ~~~
# * misc. code snippets used throughout the development process

# time-to-run timing snippet:
ST <- Sys.time()
# insert code to be timed here
ET <- Sys.time()
time_taken <- paste("expr_list reduction:", difftime(ET, ST, units = "secs"), "seconds")

# base R environment reset snippet:
ip <- as.data.frame(installed.packages())
ip <- ip[!(ip[, "Priority"] %in% c("base", "recommended")), ]
path.lib <- unique(ip$LibPath)
pkgs.to.remove <- ip[, 1]
sapply(pkgs.to.remove, remove.packages, lib = path.lib)

# mass env cleanup
rm(list = setdiff(ls(), c("KEEP VARS HERE")))
gc()

# figures.R code
# - 1.1 ~ Figure 1
# - 1.2 ~ Figure 2
# - 1.2.1 ~ 2a. BRCA2 ~ ajcc/nhg/subtype
p2a1 <- p_box(df_tcga, "ajcc_stage", "BRCA2",
              "TCGA: BRCA2 vs AJCC stage", test = "kruskal.test")
p2a2 <- p_box(df_tcga, "Subtype", "BRCA2",
              "TCGA: BRCA2 vs subtype", test = "kruskal.test")
p2a3 <- p_box(df_gse, "nhg", "BRCA2",
              "GSE: BRCA2 vs NHG", test = "kruskal.test")
p2a4 <- p_box(df_gse, "Subtype", "BRCA2",
              "GSE: BRCA2 vs subtype", test = "kruskal.test")
# - 1.2.2 ~ 2b. MKI67
p2b1 <- p_box(df_tcga, "BRCA2_group", "MKI67",
              "TCGA: MKI67 by BRCA2 group", test = "wilcox.test")
p2b2 <- p_box(df_gse, "BRCA2_group", "MKI67",
              "GSE: MKI67 by BRCA2 group", test = "wilcox.test")
# - 1.2.3 ~ stacking

# generalized p_box fn
p_box <- function(df, x, y, title = NULL, test = c("wilcox.test", "kruskal.test")) {
  test <- match.arg(test)
  ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_boxplot() +
    stat_compare_means(method = test, label.y.npc = "top") +
    theme_bw() +
    labs(title = title, x = NULL, y = NULL)
}
