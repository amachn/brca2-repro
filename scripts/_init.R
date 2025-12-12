# ~~~ _init.R ~~~
# * ONLY RUN THIS SCRIPT TO GENERATE A NEW RENV (not recommended)
# - 1 ~ renv (virtual environment)
setwd("D:/repos/brca2-repro")

install.packages("renv")
library(renv)
renv::init()

# - 2 ~ CRAN
install.packages(c(
  "msigdbr",
  "tidyverse"
))

# - 3 ~ Bioconductor
install.packages("BiocManager")
BiocManager::install(c(
  "GEOquery",
  "GSVA",
  "SummarizedExperiment",
  "org.Hs.eg.db"
))

# - 4 ~ Github
install.packages("devtools")
devtools::install_github("12379Monty/GSE96058")

# - 5 ~ Snapshots
renv::snapshot()
writeLines(capture.output(sessionInfo()), "session_info.txt")
