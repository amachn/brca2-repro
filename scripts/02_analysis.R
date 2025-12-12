# ~~~ 02_analysis.R ~~~
# - 0 ~ Libraries + Data
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(tidyverse)
})

load("data/processed/SE_tcga_gse.RData")

# --- GDC (TCGA-BRCA) ---
# - 1 ~ Extractions
# - 1.1 ~ Helper Functions
source("scripts/_utils.R")

# - 1.2 ~ Extract: Expression Values
tcga_mat <- assay(se_tcga, "log2TPM")
gse_mat <- assay(se_gse, "log2FPKM")

genes_tcga <- c(
  BRCA2="ENSG00000139618", MKI67="", E2F1="",
  RB1="", PALB2="", PARP="", ESR1=""
)

# for (g in names(genes_tcga)) colData(se_tcga)[[g]] <- tcga_mat[genes_tcga[g],]
# - for loop again for GSE, use symbols directly
colData(se_tcga)$BRCA2 <- tcga_mat[genes_tcga["BRCA2"], ]
colData(se_gse)$BRCA2 <- gse_mat["BRCA2", ]

# - 1.3 ~ Extract: AJCC + lymph nodes
extract_diag <- function(diagnoses) {
  diag <- pick_primary(diagnoses)
  patho_d <- if (!is.null(diag)) pick_first_df(diag$pathology_details) else NULL
  list(
    ajcc_stage = get_v1(diag, "ajcc_pathologic_stage"),
    ajcc_t     = get_v1(diag, "ajcc_pathologic_t"),
    ajcc_n     = get_v1(diag, "ajcc_pathologic_n"),
    ajcc_m     = get_v1(diag, "ajcc_pathologic_m"),
    ajcc_ed    = get_v1(diag, "ajcc_staging_system_edition"),
    ln_pos     = get_v1(patho_d, "lymph_nodes_positive", as_int = TRUE),
    ln_tested  = get_v1(patho_d, "lymph_nodes_tested", as_int = TRUE)
  )
}

cd <- as.data.frame(colData(se_tcga)) |>
  mutate(diag_fields = map(diagnoses, extract_diag)) |>
  unnest_wider(diag_fields)

# - 1.4 ~ Extract: ER/PR/HER2 results + expressions
cd <- cd |>
  mutate(
    ER_status   = map_chr(follow_ups, ~ get_gene_test(.x, "ESR1")),
    PR_status   = map_chr(follow_ups, ~ get_gene_test(.x, "PGR")),
    HER2_status = map_chr(follow_ups, ~ get_gene_test(.x, "ERBB2"))
  )

# - 1.x ~ extract (or determine??) mutation/HRD/heterogeneity
# note: "calculated per previously reported study by Thorsson et al."
#  - https://pubmed.ncbi.nlm.nih.gov/29628290/

# - 1.5 ~ Merge with Summarized Experiment?
# - 2 ~ Determinations
# FOR:
# - split BRCA groups? high/low
# - case_when for clinical subtype
# - hallmark pathways GSVA analysis (msigdbr, E2F, G2M, etc)
# ALSO TODO IN MORNING:
# - include master in 01_build save
# - add "ethnicity" from clinical demographics
