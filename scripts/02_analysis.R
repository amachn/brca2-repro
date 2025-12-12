# ~~~ 02_analysis.R ~~~
# - 0 ~ Libraries + Data
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(GSVA)
  library(SummarizedExperiment)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(tidyverse)
})

load("data/processed/SE_tcga_gse.RData")

# - 1 ~ Extractions
# - 1.1 ~ Helper Functions
source("scripts/_utils.R")

# - 1.2 ~ Extract: Expression Values
tcga_mat <- assay(se_tcga, "log2TPM")
gse_mat <- assay(se_gse, "log2FPKM")

genes_tcga <- c(
  BRCA2="ENSG00000139618", MKI67="ENSG00000148773", E2F1="ENSG00000101412",
  RB1="ENSG00000139687", PARP1="ENSG00000143799", # PALB2="", # doesn't exist?
  ESR1="ENSG00000091831", PGR="ENSG00000082175", ERBB2="ENSG00000141736"
)

for (g in names(genes_tcga)) {
  colData(se_tcga)[[g]] <- tcga_mat[genes_tcga[g], ]
}

for (g in c("BRCA2","MKI67","E2F1","RB1","PARP1","ESR1","PGR","ERBB2")) {
  colData(se_gse)[[g]] <- gse_mat[g, ]
}

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

# - 1.4 ~ Extract: ER/PR/HER2 results
cd <- cd |>
  mutate(
    ER_status   = map_chr(follow_ups, ~ get_gene_test(.x, "ESR1")),
    PR_status   = map_chr(follow_ups, ~ get_gene_test(.x, "PGR")),
    HER2_status = map_chr(follow_ups, ~ get_gene_test(.x, "ERBB2"))
  )

# - 1.x ~ extract (or determine??) mutation/HRD/heterogeneity
# note: "calculated per previously reported study by Thorsson et al."
#  - https://pubmed.ncbi.nlm.nih.gov/29628290/

# - 1.5 ~ Merge with Summarized Experiment
colData(se_tcga) <- DataFrame(cd)

# - 2 ~ Determinations
# - 2.1 ~ Calculate BRCA Groups
colData(se_tcga)$BRCA2_group <- ifelse(
  colData(se_tcga)$BRCA2 > median(colData(se_tcga)$BRCA2, na.rm = TRUE),
  "High", "Low"
)
colData(se_gse)$BRCA2_group <- ifelse(
  colData(se_gse)$BRCA2 > median(colData(se_gse)$BRCA2, na.rm = TRUE),
  "High", "Low"
)

# - 2.2 ~ Hallmark Pathway Scores [msigdbr, E2F, G2M, etc] (GSVA)
# - 2.2.1 ~ TCGA Mapping (Ensembl IDs -> Gene Symbols)
ens_ids <- rownames(tcga_mat)
symbols <- mapIds(
  org.Hs.eg.db, keys = ens_ids, column = "SYMBOL",
  keytype = "ENSEMBL", multiVals = "first"
)

keep <- !is.na(symbols) # drop unmapped genes
tcga_mat_sym <- tcga_mat[keep, ]
rownames(tcga_mat_sym) <- symbols[keep]
tcga_mat_sym <- rowsum(tcga_mat_sym, group = rownames(tcga_mat_sym)) # collapse

# - 2.2.2 ~ Run GSVA
msig <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- split(msig$gene_symbol, msig$gs_name)

tcga_scores <- gsva(gsvaParam(
  exprData = tcga_mat_sym, geneSets = hallmark, kcdf = "Gaussian"
))
gse_scores <- gsva(gsvaParam(
  exprData = gse_mat, geneSets = hallmark, kcdf = "Gaussian"
))

# - 2.2.3 ~ Store Pathway Scores
colData(se_tcga)$E2F <- tcga_scores["HALLMARK_E2F_TARGETS", ]
colData(se_tcga)$G2M <- tcga_scores["HALLMARK_G2M_CHECKPOINT", ]
colData(se_tcga)$P53 <- tcga_scores["HALLMARK_P53_PATHWAY", ]
colData(se_tcga)$ER_early <- tcga_scores["HALLMARK_ESTROGEN_RESPONSE_EARLY", ]
colData(se_tcga)$ER_late <- tcga_scores["HALLMARK_ESTROGEN_RESPONSE_LATE", ]

colData(se_gse)$E2F <- gse_scores["HALLMARK_E2F_TARGETS", ]
colData(se_gse)$G2M <- gse_scores["HALLMARK_G2M_CHECKPOINT", ]
colData(se_gse)$P53 <- gse_scores["HALLMARK_P53_PATHWAY", ]
colData(se_gse)$ER_early <- gse_scores["HALLMARK_ESTROGEN_RESPONSE_EARLY", ]
colData(se_gse)$ER_late <- gse_scores["HALLMARK_ESTROGEN_RESPONSE_LATE", ]

# - 2.2.4 ~ sanity check
summary(colData(se_tcga)$P53)
summary(colData(se_gse)$P53)

# - 2.3 ~ Clinical Subtype Determination (case_when)
# - 2.3.1 ~ TCGA
cd_tcga <- as.data.frame(colData(se_tcga))
cd_tcga$Subtype <- case_when(
  cd_tcga$ER_status == "Positive" ~ "ER+",
  cd_tcga$HER2_status == "Positive" ~ "HER2+",
  cd_tcga$ER_status == "Negative" &
    cd_tcga$PR_status == "Negative" &
    cd_tcga$HER2_status == "Negative" ~ "TNBC",
  TRUE ~ NA_character_
)
colData(se_tcga)$Subtype <- cd_tcga$Subtype

# - 2.3.2 ~ GSE
cd_gse <- as.data.frame(colData(se_gse))
cd_gse$Subtype <- case_when(
  cd_gse$er_Status == "1" ~ "ER+",
  cd_gse$her2_Status == "1" ~ "HER2+",
  cd_gse$er_Status == "0" &
    cd_gse$pgr_Status == "0" &
    cd_gse$her2_Status == "0" ~ "TNBC",
  TRUE ~ NA_character_
)
colData(se_gse)$Subtype <- cd_gse$Subtype

# - 2.3.3 ~ sanity check
table(colData(se_tcga)$Subtype)
table(colData(se_gse)$Subtype)

# - ext ~ env cleanup
rm(list = setdiff(
  ls(), c("se_tcga", "se_gse", "master", "tcga_scores", "gse_scores")
))
gc()

# ~~~ OUTPUT -> 03_figures.R ~~~
save(se_tcga, se_gse, master, tcga_scores, gse_scores,
     file = "data/processed/SE_tcga_gse_anno.RData")
