# Check to make sure PharmacoGx is installed

if (!requireNamespace("PharmacoGx", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("PharmacoGx")
}


# Load required library
library(PharmacoGx)

# ---------- CONFIGURATION ----------
datasets <- c(
  "FIMM_2016", "BeatAML_2018", "CCLE_2015", "CTRPv2_2015", "gCSI_2019",
  "GDSC_2020(v1-8.2)", "GDSC_2020(v2-8.2)", "NCI60_2021", "PDTX_2019", "PRISM_2020"
)

pset_dir <- "PharmacoGx/pharmacogx_psets"
export_dir <- "PharmacoSet_Exports"

dir.create(pset_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- HELPER FUNCTION: Save DataFrame with Row Names as a Column ----------
write_with_rownames <- function(df, file, rowname_col = "ID") {
  df_out <- cbind(setNames(data.frame(rownames(df)), rowname_col), df)
  write.csv(df_out, file, row.names = FALSE)
}

# ---------- HELPER FUNCTION: Flatten Raw Sensitivity Data ----------
export_sensitivity_raw_flat <- function(pset, out_file, id_col = "ExperimentID") {
  raw_array <- sensitivityRaw(pset)
  info_df <- sensitivityInfo(pset)

  rows <- vector("list", length = dim(raw_array)[1])

  for (i in seq_len(dim(raw_array)[1])) {
    dose_vals <- raw_array[i, , "Dose"]
    viability_vals <- raw_array[i, , "Viability"]

    max_doses <- length(dose_vals)
    dose_names <- paste0("dose.", seq_len(max_doses), ".Dose")
    viability_names <- paste0("dose.", seq_len(max_doses), ".Viability")

    row_data <- setNames(as.list(c(dose_vals, viability_vals)), c(dose_names, viability_names))
    meta_info <- as.list(info_df[i, ])
    rows[[i]] <- c(meta_info, row_data)
  }

  final_df <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  rownames(final_df) <- rownames(info_df)

  final_df <- cbind(setNames(data.frame(rownames(final_df)), id_col), final_df)
  write.csv(final_df, out_file, row.names = FALSE)
}

# ---------- MAIN LOOP ----------
for (dataset in datasets) {
  cat("\nProcessing:", dataset, "\n")

  # Construct paths
  pset_path <- file.path(pset_dir, paste0(dataset, ".rds"))
  out_dir <- file.path(export_dir, dataset)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ---------- DOWNLOAD STEP ----------
  if (!file.exists(pset_path)) {
    cat("Downloading PSet from PharmacoDB:", dataset, "\n")
    tryCatch({
      pset <- downloadPSet(dataset)
      saveRDS(pset, pset_path)
      cat("Downloaded and saved:", pset_path, "\n")
    }, error = function(e) {
      cat("Failed to download:", dataset, "->", e$message, "\n")
      next
    })
  } else {
    pset <- tryCatch(readRDS(pset_path), error = function(e) {
      cat("Failed to read:", dataset, "->", e$message, "\n")
      return(NULL)
    })
  }
  if (is.null(pset)) next

  # ---------- 1. SENSITIVITY DATA ----------

  try({
    export_sensitivity_raw_flat(pset, file.path(out_dir, paste0(dataset, "_sensitivity_raw.csv")))
    cat("Saved flattened sensitivity raw data\n")
  }, silent = TRUE)

  try({
    summary <- summarizeSensitivityProfiles(pset, summary.stat = "mean")
    write_with_rownames(summary, file.path(out_dir, paste0(dataset, "_sensitivity_summary.csv")), "ExperimentID")
    cat("Saved sensitivity summary\n")
  }, silent = TRUE)

  try({
    info <- sensitivityInfo(pset)
    write_with_rownames(info, file.path(out_dir, paste0(dataset, "_sensitivity_info.csv")), "ExperimentID")
    cat("Saved sensitivity info\n")
  }, silent = TRUE)

  try({
    profiles <- sensitivityProfiles(pset)
    write_with_rownames(profiles, file.path(out_dir, paste0(dataset, "_sensitivity_profiles.csv")), "ExperimentID")
    cat("Saved sensitivity profiles\n")
  }, silent = TRUE)

  # ---------- 2. CELL & DRUG METADATA ----------

  try({
    write_with_rownames(cellInfo(pset), file.path(out_dir, paste0(dataset, "_cellinfo.csv")), "CellLine")
    cat("Saved cell info\n")
  }, silent = TRUE)

  try({
    write_with_rownames(drugInfo(pset), file.path(out_dir, paste0(dataset, "_druginfo.csv")), "Drug")
    cat("Saved drug info\n")
  }, silent = TRUE)

  # ---------- 3. MOLECULAR (OMICS) DATA ----------

  omics_types <- names(molecularProfiles(pset))
  cat("Omics types detected:", paste(omics_types, collapse = ", "), "\n")

  for (datatype in omics_types) {
    cat("Exporting:", datatype, "\n")

    mat <- tryCatch({
      data_obj <- molecularProfiles(pset)[[datatype]]
      if (inherits(data_obj, "SummarizedExperiment")) {
        assay(data_obj)
      } else if (is.matrix(data_obj) || is.data.frame(data_obj)) {
        data_obj
      } else {
        stop("Unsupported omics data type")
      }
    }, error = function(e) {
      cat("Error extracting", datatype, ":", conditionMessage(e), "\n")
      return(NULL)
    })

    if (!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0) {
      outfile <- file.path(out_dir, paste0(dataset, "_molecular_", datatype, ".csv"))
      write_with_rownames(mat, outfile, rowname_col = "Feature")
      cat("Saved", datatype, "(", nrow(mat), "×", ncol(mat), ")\n")

      # binary mutation export
      if (datatype %in% c("mutation", "mutationall", "mutationchosen")) {
        if (nrow(mat) * ncol(mat) < 1e7) {
          mut_bin <- as.data.frame(mat != "wt")
          mut_bin[] <- lapply(mut_bin, as.integer)
          write_with_rownames(mut_bin, file.path(out_dir, paste0(dataset, "_molecular_", datatype, "_binary.csv")), "Feature")
          cat("Saved binary mutation matrix\n")
        } else {
          cat("Skipped binary mutation (too large)\n")
        }
      }
    } else {
      cat("Skipped", datatype, "- matrix empty or invalid\n")
    }

    rm(mat)
    gc()
  }

  # Clean up memory
  rm(pset)
  gc()
}

cat("\n All datasets processed.\n")
