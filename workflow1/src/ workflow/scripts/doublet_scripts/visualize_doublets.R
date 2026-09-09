# scripts/doublet_scripts/visualize_doublets_one_sample.R

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(future)
})

input_file      <- snakemake@input[["seurat_with_doublets"]]
output_plot     <- snakemake@output[["umap_plot_single"]]
log_file        <- snakemake@log[[1]]
confounders     <- snakemake@params[["confounders"]]
n_pcs           <- snakemake@params[["n_pcs"]]
sample_id       <- snakemake@wildcards[["sample_id"]]

if (!dir.exists(dirname(log_file))) dir.create(dirname(log_file), recursive = TRUE)
log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] [%s] %s", timestamp, level, sample_id, msg)
  cat(line, "\n", file = log_file, append = TRUE)
}

tryCatch({
  log_msg("Script starts execution (visualized two cells).")
  plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)
  
  log_msg(paste("Loading Seurat object with doublet labels:", basename(input_file)))
  seurat_obj <- readRDS(input_file)
  log_msg(paste("Object contains", ncol(seurat_obj), "cells."))
  
  log_msg("Performing temporary normalization for visualization...")
  seurat_obj <- tryCatch({
    SCTransform(
      seurat_obj,
      assay = "RNA",
      vars.to.regress = confounders,
      verbose = FALSE,
      conserve.memory = TRUE
    )
  }, error = function(e) {
    log_msg(paste("SCTransform failed, falling back to NormalizeData workflow:", e$message), level = "WARN")
    x <- NormalizeData(seurat_obj, verbose = FALSE)
    x <- FindVariableFeatures(x, nfeatures = 2000, verbose = FALSE)
    x <- ScaleData(x, vars.to.regress = confounders, verbose = FALSE)
    x
  })
  
  log_msg("Running PCA and UMAP...")
  seurat_obj <- RunPCA(seurat_obj, npcs = n_pcs, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs, verbose = FALSE)
  
  doublet_rate <- mean(seurat_obj$doublet_class == "doublet") * 100
  
  log_msg("Generating UMAP plot...")
  color_map <- c("singlet" = "grey", "doublet" = "red")
  
  umap_plot <- DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "doublet_class",
    cols = color_map,
    pt.size = 0.5
  ) +
    ggtitle(paste("Sample:", sample_id), 
            subtitle = sprintf("Doublet Rate: %.2f%%", doublet_rate)) +
    theme_bw()
  
  ggsave(output_plot, plot = umap_plot, width = 7, height = 6)
  log_msg(paste("Single sample UMAP plot saved to:", basename(output_plot)))
  
  log_msg("Script executed successfully.")
  
}, error = function(e) {
  log_msg(paste("ERROR:", e$message), level = "ERROR")
  stop(e)
})
