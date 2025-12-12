stopifnot(
  requireNamespace("dplyr", quietly = TRUE),
  requireNamespace("purrr", quietly = TRUE)
)

# safely get first non-null dataframe, or return NULL
pick_first_df <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x)) {
    x <- compact(x)
    if (length(x)) return(as.data.frame(x[[1]]))
  }
  NULL
}

# safely get first value of a column
get_v1 <- function(df, col, as_int = FALSE) {
  if (is.null(df) || nrow(df) == 0 || !col %in% names(df))
    return(if (as_int) NA_integer_ else NA_character_)
  out <- df[[col]][1]
  if (as_int) suppressWarnings(as.integer(out)) else as.character(out)
}

# select the relevant primary row (nested dataframes); falls back if needed
pick_primary <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  i <- if ("classification_of_tumor" %in% names(df)) {
    match("primary", df$classification_of_tumor)
  } else NA_integer_
  if (!is.na(i)) df[i, , drop = FALSE] else df[1, , drop = FALSE]
}

# return mol. test result for a gene symbol from follow_ups
get_gene_test <- function(fu_df, gene) {
  if (is.null(fu_df) || nrow(fu_df) == 0 || !"molecular_tests" %in% names(fu_df))
    return(NA_character_)

  tests_df <- compact(fu_df$molecular_tests)
  if (!length(tests_df)) return(NA_character_)
  tests_df <- bind_rows(tests_df)

  if (!all(c("gene_symbol", "test_result") %in% names(tests_df)))
    return(NA_character_)

  vals <- tests_df |>
    filter(gene_symbol == gene, !is.na(test_result)) |>
    pull(test_result)

  if (!length(vals)) return(NA_character_)
  if (any(vals == "Positive")) return("Positive")
  if (any(vals == "Equivocal")) return("Equivocal")
  if (any(vals == "Negative")) return("Negative")
  sub[1]
}
