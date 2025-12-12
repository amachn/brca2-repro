# ~~~ 02_analysis.R ~~~
# - 0 ~ Libraries + Data
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
})

load("data/processed/SE_tcga_gse.RData")