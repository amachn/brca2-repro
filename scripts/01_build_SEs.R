# ~~~ 01_build_SEs.R ~~~
# - 0 ~ Libraries
suppressPackageStartupMessages({
  # library(GEOquery) # only needed if you plan on running getGEOSuppFiles
  library(SummarizedExperiment)
  library(dplyr)
  library(jsonlite)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

# --- GDC (TCGA-BRCA) ---
# - 1 ~ Folder/File Setup
setwd("D:/repos/brca2-repro")

data_dir <- "data/raw/gdc_data"
expr_dir <- file.path(data_dir, "expression")

metadata_json <- file.path(data_dir, "metadata.json")
clinical_json <- file.path(data_dir, "clinical.json")
sample_tsv <- file.path(data_dir, "sample_sheet.tsv")

# - 2 ~ Read Files
meta <- fromJSON(metadata_json, flatten = TRUE) |>
  as_tibble() |>
  unnest(associated_entities) |>
  transmute(
    file_id = file_id,
    file_name = file_name,
    case_id = case_id,
    ent_id = entity_id,
    entSub_id = entity_submitter_id
  )

clinical <- fromJSON(clinical_json, flatten = TRUE) |>
  as_tibble() |>
  transmute(
    case_id = case_id,
    submitter_id = submitter_id,
    disease_type = disease_type,
    age_at_diag = demographic.age_at_index,
    vital_status = demographic.vital_status,
    ethnicity = demographic.ethnicity,
    race = demographic.race,
    diagnoses = diagnoses,
    follow_ups = follow_ups
  )

samp <- read_tsv(sample_tsv, show_col_types = FALSE) |>
  rename_with(~ str_replace_all(., " ", "_")) |>
  rename_with(str_to_lower) |>
  dplyr::rename(
    samp_case_id = case_id
  )

# - 3 ~ Join into Master Table
master <- samp |>
  left_join(meta, by = c("file_id", "file_name")) |>
  left_join(clinical, by = "case_id")

# - 4 ~ Read all Expression Files
read_expr_file <- function(path) {
  read_tsv(path, comment = "#", show_col_types = FALSE) |>
    filter(!grepl("^N_", gene_id)) |>
    select(gene_id, tpm_unstranded)
}

identifiers <- master |>
  select(file_id, file_name, sample_id) |>
  distinct()

# may take a while to run (avg. 1-2 min)
expr_list <- purrr::map(seq_len(nrow(identifiers)), function(i) {
  folder_name <- identifiers$file_id[i]
  file_name <- identifiers$file_name[i]
  sid <- identifiers$sample_id[i]
  fp <- file.path(expr_dir, folder_name, file_name)

  if (!file.exists(fp)) {
    warning("Missing expr file: ", fp)
    return(NULL)
  }

  read_expr_file(fp) |>
    dplyr::rename(!!sid := tpm_unstranded)
})

# may take a while to run (avg. 2-5 min)
expr_tbl <- purrr::reduce(expr_list, full_join, by = "gene_id")

expr_mat <- expr_tbl |>
  mutate(gene_id = str_replace(gene_id, "\\.\\d+$", "")) |>
  distinct(gene_id, .keep_all = TRUE) |>
  column_to_rownames("gene_id") |>
  as.matrix()

expr_tcga_log2 <- log2(expr_mat + 0.1) # match GSE96058 scale

# - 5 ~ Align TCGA Expr w/ Master + Summarize Experiment
common_samples <- intersect(colnames(expr_tcga_log2), master$sample_id)
expr_mat_sub <- expr_tcga_log2[, common_samples, drop = FALSE]
sample_data <- master |>
  filter(sample_id %in% common_samples) |>
  arrange(match(sample_id, colnames(expr_mat_sub)))

se_tcga <- SummarizedExperiment(
  assays = list(log2TPM = expr_mat_sub),
  colData = as.data.frame(sample_data)
)

# - extra ~ GDC Memory Cleanup
rm(list = c(
  "clinical_json", "data_dir", "expr_dir", "metadata_json", "sample_tsv",
  "meta", "clinical", "samp", "read_expr_file", "identifiers", # "master",
  "expr_list", "expr_tbl", "expr_mat", "expr_tcga_log2", "common_samples",
  "expr_mat_sub", "sample_data"
))
gc()

# --- GEO (GSE96058) ---
# - 1 ~ Download -> Load GSE Expression Data
# getGEOSuppFiles("GSE96058", makeDirectory = TRUE) # one-time run!
raw_gse_expr <- read_csv( # renamed & extracted from ^
  "data/raw/gse_data/GSE_GEXP_3273s_136rt.csv", show_col_types = FALSE)

gse_expr <- as.matrix(raw_gse_expr[, -1])
rownames(gse_expr) <- raw_gse_expr[[1]]

# - 2 ~ Load GSE Clinical Data
data("GSE96058_sampDesc", package = "GSE96058")

sampDesc <- as.data.frame(GSE96058_sampDesc) |> filter(!isRepl)

# - 3 ~ Align GSE Expr w/ Clinical Data + Summarize Experiment
meta_ids <- sampDesc$title
common_ids <- intersect(colnames(gse_expr), meta_ids)
samp_rows <- match(common_ids, meta_ids)

gse_expr_aligned <- gse_expr[, common_ids, drop = FALSE]
sampDesc_aligned <- sampDesc[samp_rows, , drop = FALSE]
rownames(sampDesc_aligned) <- common_ids

se_gse <- SummarizedExperiment(
  assays = list(log2FPKM = gse_expr_aligned),
  colData = sampDesc_aligned
)

# - ext ~ GEO Cleanup
rm(list = c(
  "raw_gse_expr", "gse_expr", "GSE96058_sampDesc", "sampDesc", "meta_ids",
  "common_ids", "samp_rows", "gse_expr_aligned", "sampDesc_aligned"
))
gc()

# ~~~ OUTPUT -> 02_analysis.R ~~~
save(se_tcga, se_gse, master, file = "data/processed/SE_tcga_gse.RData")
