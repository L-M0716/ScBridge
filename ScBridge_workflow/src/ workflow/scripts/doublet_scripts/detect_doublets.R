# scripts/doublet_scripts/detect_doublets_one_sample.R

suppressPackageStartupMessages({
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
})

input_file      <- snakemake@input[["seurat_qc"]]
output_seurat   <- snakemake@output[["seurat_with_doublets"]]
output_stats    <- snakemake@output[["doublet_stats"]]
expected_rate   <- snakemake@params[["expected_rate"]]
seed            <- snakemake@params[["seed"]]
method          <- snakemake@params[["method"]]
log_file        <- snakemake@log[[1]]
sample_id       <- snakemake@wildcards[["sample_id"]]

if (!dir.exists(dirname(log_file))) dir.create(dirname(log_file), recursive = TRUE)
log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] [%s] %s", timestamp, level, sample_id, msg)
  cat(line, "\n", file = log_file, append = TRUE)
}

tryCatch({
  log_msg("Script started (Doublet Detection).")
  log_msg(paste("Loading QC'd Seurat object:", basename(input_file)))
  seurat_obj <- readRDS(input_file)
  log_msg(paste("Object contains", ncol(seurat_obj), "cells."))
  

  log_msg("Converting Seurat object to SingleCellExperiment object...")
  sce <- as.SingleCellExperiment(seurat_obj)

  set.seed(seed)
  log_msg(paste("Using", method, "for detection, expected doublet rate:", expected_rate))
  sce <- scDblFinder(sce, dbr = expected_rate)
  

  log_msg("Converting doublet detection results back to Seurat object...")
  seurat_obj$doublet_score <- sce$scDblFinder.score
  seurat_obj$doublet_class <- sce$scDblFinder.class

  seurat_obj$doublet_call  <- ifelse(sce$scDblFinder.class == "doublet", 1, 0)
  
  log_msg("Generating statistics...")
  stats_df <- data.frame(
    sample_id = sample_id,
    total_cells = ncol(seurat_obj),
    doublet_count = sum(seurat_obj$doublet_call),
    doublet_rate_percent = mean(seurat_obj$doublet_call) * 100,
    detection_method = method
  )
  
  if (!dir.exists(dirname(output_seurat))) dir.create(dirname(output_seurat), recursive = TRUE)
  saveRDS(seurat_obj, file = output_seurat)
  log_msg(paste("QC'd Seurat object with doublet annotations saved to:", basename(output_seurat)))
  
  if (!dir.exists(dirname(output_stats))) dir.create(dirname(output_stats), recursive = TRUE)
  write.csv(stats_df, file = output_stats, row.names = FALSE)
  log_msg(paste("Doublet statistics saved to:", basename(output_stats)))
  log_msg("Script executed successfully.")
  
}, error = function(e) {
  log_msg(paste("ERROR:", e$message), level = "ERROR")
  stop(e)
})