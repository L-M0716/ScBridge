library(Seurat)
library(dplyr)

seurat_input     <- snakemake@input[["seurat_raw"]]
qc_table_output  <- snakemake@output[["qc_table"]]
seurat_qc_output <- snakemake@output[["seurat_qc"]]
min_genes        <- snakemake@params[["min_genes"]]
max_genes        <- snakemake@params[["max_genes"]]
min_counts       <- snakemake@params[["min_counts"]]
max_mito         <- snakemake@params[["max_mito"]]
log_file         <- snakemake@log[[1]]
sample_id        <- snakemake@wildcards[["sample_id"]]

if (!dir.exists(dirname(log_file))) {
  dir.create(dirname(log_file), recursive = TRUE)
}
log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] [%s] %s", timestamp, level, sample_id, msg)
  cat(line, "\n", file = log_file, append = TRUE)
}

tryCatch({
  log_msg("Script execution started (calculating QC metrics)")
  

  seurat_obj <- readRDS(seurat_input)
  

  seurat_obj$qc_status <- ifelse(
    seurat_obj$nFeature_RNA >= min_genes &
      seurat_obj$nFeature_RNA <= max_genes &
      seurat_obj$nCount_RNA >= min_counts &
      seurat_obj$percent.mt <= max_mito,
    "pass", "fail"
  )
  

  if (!dir.exists(dirname(qc_table_output))) dir.create(dirname(qc_table_output), recursive = TRUE)
  write.csv(seurat_obj@meta.data, file = qc_table_output, row.names = TRUE)
  log_msg(paste("The QC indicator table has been saved to:", basename(qc_table_output)))
  

  saveRDS(seurat_obj, file = seurat_qc_output)
  log_msg(paste("The Seurat object with QC markers has been saved to:", basename(seurat_qc_output)))
  
  log_msg("Script execution completed successfully.")
  
}, error = function(e) {
  log_msg(paste("ERROR:", e$message), level = "ERROR")
  stop(e$message)
})