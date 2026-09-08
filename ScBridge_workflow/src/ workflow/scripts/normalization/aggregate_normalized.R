# scripts/normalization_scripts/aggregate_normalized.R
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

normalized_files <- snakemake@input[["normalized_list"]]
output_merged    <- snakemake@output[["merged_normalized_seurat"]]
output_qc_fig    <- snakemake@output[["qc_figures"]]
output_hvg       <- snakemake@output[["hvg_list"]]
n_variable_genes <- snakemake@params[["n_variable_genes"]]
log_file         <- snakemake@log[[1]]

if (!dir.exists(dirname(log_file))) dir.create(dirname(log_file), recursive = TRUE)
log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] [%s] %s", timestamp, level, "aggregate", msg)
  cat(line, "\n", file = log_file, append = TRUE)
}

tryCatch({
  log_msg("Script started (aggregating normalized results).")
  
  log_msg(paste("Reading", length(normalized_files), "normalized Seurat objects..."))
  seurat_list <- lapply(normalized_files, readRDS)
  seurat_list <- lapply(seurat_list, function(x) {
    tryCatch({
      DietSeurat(
        x,
        assays = c("RNA", "SCT"),
        dimreducs = NULL,
        graphs = NULL,
        misc = FALSE,
        layers = "counts"
      )
    }, error = function(e) {
      x
    })
  })
  

  log_msg(sprintf("Selecting %d consensus highly variable genes from all samples...", n_variable_genes))
  integration_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = n_variable_genes)
  

  log_msg("Starting to merge all normalized Seurat objects...")
  sample_ids <- gsub("/seurat_normalized.rds", "", basename(dirname(normalized_files)))
  names(seurat_list) <- sample_ids
  
  if (length(seurat_list) == 1) {
    merged_seurat <- seurat_list[[1]]
  } else {
    merged_seurat <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)], add.cell.ids = names(seurat_list))
  }
  

  DefaultAssay(merged_seurat) <- "SCT"
  VariableFeatures(merged_seurat) <- integration_features
  log_msg(sprintf("Merge completed, final object contains %d cells. Default assay set to SCT, highly variable genes set to %d.",
                  ncol(merged_seurat), length(integration_features)))
  

  p1 <- tryCatch({
    VlnPlot(merged_seurat, features = "nCount_SCT", group.by = "sample_id", pt.size = 0) + 
      labs(title = "Normalized UMI Distribution by Sample")
  }, error = function(e) {
    log_msg(paste("VlnPlot skipped:", e$message), level = "WARN")
    NULL
  })
  

  if (!dir.exists(dirname(output_merged))) dir.create(dirname(output_merged), recursive = TRUE)
  saveRDS(merged_seurat, file = output_merged, compress = FALSE)
  log_msg(paste("Final merged object saved to:", basename(output_merged)))
  
  writeLines(integration_features, con = output_hvg)
  log_msg(paste("Consensus highly variable genes list saved to:", basename(output_hvg)))
  
  if (!is.null(p1)) {
    ggsave(output_qc_fig, plot = p1, width = 8, height = 6)
    log_msg(paste("Normalized QC plot saved to:", basename(output_qc_fig)))
  } else {
    pdf(output_qc_fig, width=8, height=6)
    plot.new()
    text(0.5, 0.5, "VlnPlot unavailable (Seurat 5 compatibility)")
    dev.off()
    log_msg("Normalized QC plot saved to placeholder.", level = "WARN")
  }
  
  log_msg("Script executed successfully.")
}, error = function(e) {
  log_msg(paste("ERROR:", e$message), level = "ERROR")
  stop(e)
})
