library(PharmacoGx)

# ---- Paths ----
input_dir  <- "pharmacogx_psets"     # folder with .rds files
output_dir <- "PharmacoSet_Exports/gene_expr"  # folder to save CSVs
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Explicitly list the datasets you want ----
datasets <- c(
  "CCLE_2015.rds",
  "gCSI_2019.rds", "GDSC_2020(v1-8.2).rds", "GDSC_2020(v2-8.2).rds"
)

# ---- Loop over datasets ----
for (dataset in datasets) {
  pset_path <- file.path(input_dir, dataset)
  fname <- basename(pset_path)
  dataset_name <- tools::file_path_sans_ext(fname)
  output_file <- file.path(output_dir, paste0(dataset_name, "_rnaseq_summarized.csv"))

  cat("\nLoading", dataset_name, "dataset...\n")
  pset <- readRDS(pset_path)

  # ---- Check molecular data types ----
  omics_types <- names(molecularProfiles(pset))
  cat("Available omics types:", paste(omics_types, collapse = ", "), "\n")

  # ---- Require Kallisto rnaseq ----
  mDataType <- "Kallisto_0.46.1.rnaseq"
  if (!(mDataType %in% omics_types)) {
    cat("Skipping — no", mDataType, "found.\n")
    next
  }

  se <- molecularProfiles(pset)[[mDataType]]

  # ---- Validate data exists ----
  if (!inherits(se, "SummarizedExperiment") || nrow(se) == 0 || ncol(se) == 0) {
    cat("Skipping —", mDataType, "assay is empty.\n")
    next
  }

  # ---- Extract sample IDs ----
  sample_ids <- colData(se)$sampleid
  if (is.null(sample_ids) || all(is.na(sample_ids))) {
    cat("Skipping — no sampleid metadata found.\n")
    next
  }

  valid_cells <- intersect(sample_ids, cellNames(pset))
  features <- rownames(se)

  if (length(valid_cells) == 0) {
    cat("Skipping — no matching cell lines found.\n")
    next
  }

  # ---- Summarize molecular profiles ----
  cat("Summarizing", mDataType, "data for", length(valid_cells), "cell lines...\n")
  summarized_expr <- summarizeMolecularProfiles(
    object       = pset,
    mDataType    = mDataType,
    cell.lines   = valid_cells,
    features     = features,
    summary.stat = "mean",
    fill.missing = TRUE,
    verbose      = TRUE
  )

  # ---- Save output ----
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  if (inherits(summarized_expr, "ExpressionSet")) {
    summarized_expr <- Biobase::exprs(summarized_expr)
  } else if (inherits(summarized_expr, "SummarizedExperiment")) {
    summarized_expr <- SummarizedExperiment::assay(summarized_expr)
  }

  write.csv(summarized_expr, output_file)
  cat("Done! File saved to:", output_file, "\n")
}