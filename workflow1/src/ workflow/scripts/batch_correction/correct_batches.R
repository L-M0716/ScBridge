#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(sva)
  library(BiocParallel)
})
options(future.globals.maxSize = 20 * 1024^3)
input_seurat   <- snakemake@input[["seurat_normalized"]]
output_corrected <- snakemake@output[["corrected"]]
output_umap      <- snakemake@output[["umap_plot"]]
output_pca       <- snakemake@output[["pca_plot"]]
method         <- snakemake@params[["method"]]
batch_key      <- snakemake@params[["batch_key"]]
n_pcs          <- snakemake@params[["n_pcs"]]
log_file       <- snakemake@log[[1]]
n_features_integrate <- snakemake@config[["integration"]][["n_features"]]
seed <- snakemake@params[["seed"]]


if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] %s", timestamp, level, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  if (level == "ERROR") stop(line)
}

validate_params <- function(seurat_obj, batch_key, n_features_integrate) {
  if (!batch_key %in% colnames(seurat_obj@meta.data)) {
    log_msg(sprintf("WARNING: Batch variable does not exist in the metadata '%s', skip batch normalization.", batch_key), "WARNING")
    return(list(seurat_obj = seurat_obj, skip_correction = TRUE))
  }
  
  batch_levels <- unique(seurat_obj[[batch_key, drop = TRUE]])
  if (length(batch_levels) <= 1) {
    log_msg(sprintf("WARNING: Batch variable '%s' only contains %d batches, skipping batch correction.", 
                    batch_key, length(batch_levels)), "WARNING")
    return(list(seurat_obj = seurat_obj, skip_correction = TRUE))
  }
  
  if ("SCT" %in% Assays(seurat_obj)) {
    DefaultAssay(seurat_obj) <- "SCT"
    log_msg("Detected SCT assay, will use SCT for downstream analysis.")
    if (length(VariableFeatures(seurat_obj)) == 0) {
      log_msg("No variable features found in SCT assay, re-running FindVariableFeatures.", "WARNING")
      seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = n_features_integrate)
    }
  } else {
    log_msg("No SCT assay detected, will use RNA assay. Please ensure data is normalized.", "WARNING")
    DefaultAssay(seurat_obj) <- "RNA"
    if (length(VariableFeatures(seurat_obj)) == 0) {
      log_msg("No variable features found in RNA assay, re-running NormalizeData and FindVariableFeatures.", "WARNING")
      seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
      seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = n_features_integrate, verbose = FALSE)
    }
  }
  
  return(list(seurat_obj = seurat_obj, skip_correction = FALSE))
}


run_correction <- function(seurat_obj, method, batch_key, n_pcs, n_features_integrate) {
  
  method_lower <- tolower(method)
  
  if (method_lower == "harmony") {
    log_msg("Using Harmony for batch correction...")

    seurat_obj <- RunHarmony(
      object = seurat_obj,
      group.by.vars = batch_key,
      reduction.use = "pca",
      dims.use = 1:n_pcs,
      verbose = FALSE
    )
    seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:n_pcs, reduction.name = "umap.harmony", verbose = FALSE,seed.use = seed)
    return(seurat_obj)
    
  } else if (method_lower == "seurat") {
    log_msg("Using Seurat v5 SCT integration workflow...")


    log_msg("Joining layered counts for SCT assay (compatible with Seurat v5 layer structure)...")
    if ("SCT" %in% names(seurat_obj@assays)) {
      tryCatch({
        seurat_obj[["SCT"]] <- JoinLayers(seurat_obj[["SCT"]])
        log_msg("SCT JoinLayers completed.")
      }, error = function(e) {
        log_msg(sprintf("Skipping JoinLayers (likely already single layer): %s", e$message), "WARNING")
      })
    }

    log_msg(sprintf("Splitting object by '%s' for integration...", batch_key))
    seurat_list <- SplitObject(seurat_obj, split.by = batch_key)

    integ_features <- VariableFeatures(seurat_obj)
    log_msg(sprintf("Using pre-selected %d integration features.", length(integ_features)))

    log_msg("Preparing SCT integration...")
    seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = integ_features, verbose = FALSE)
    
    log_msg("Finding integration anchors...")
    integ_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = integ_features, verbose = FALSE)
    
    log_msg("Integrating data...")
    seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT", verbose = FALSE)
    if (!"RNA" %in% names(seurat_integrated@assays) && "RNA" %in% names(seurat_obj@assays)) {
      log_msg("Detected RNA assay lost after integration, re-attaching from original object...")
      seurat_integrated[["RNA"]] <- seurat_obj[["RNA"]]
      seurat_integrated <- NormalizeData(seurat_integrated, assay = "RNA",
                                         normalization.method = "LogNormalize",
                                         scale.factor = 10000, verbose = FALSE)
      log_msg("The RNA assay has been remounted and log-normalization has been completed.")
    }
    log_msg("Performing dimensionality reduction on integrated data...")
    DefaultAssay(seurat_integrated) <- "integrated"
    seurat_integrated <- RunPCA(seurat_integrated, npcs = n_pcs, verbose = FALSE, seed.use = seed)
    seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:n_pcs, reduction.name = "umap.integrated", verbose = FALSE, seed.use = seed)
    
    return(seurat_integrated)
    
  } else if (method_lower == "combat") {
    log_msg("Using sva::ComBatSeq for batch correction...")
    DefaultAssay(seurat_obj) <- "RNA"
    counts_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
    meta_data <- seurat_obj@meta.data
    
    corrected_counts <- ComBat_seq(
      counts = as.matrix(counts_matrix), 
      batch = meta_data[[batch_key]], 
      group = NULL
    )
    
    seurat_obj[["ComBat"]] <- CreateAssayObject(counts = corrected_counts)
    DefaultAssay(seurat_obj) <- "ComBat"
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE) %>%
      FindVariableFeatures(verbose = FALSE) %>%
      ScaleData(verbose = FALSE) %>%
      RunPCA(npcs = n_pcs, verbose = FALSE) %>%
      RunUMAP(dims = 1:n_pcs, reduction.name = "umap.combat", verbose = FALSE)
    return(seurat_obj)
    
  } else {
    log_msg(sprintf("Unsupported correction method: %s", method), "ERROR")
  }
}

generate_plots <- function(obj, reduction_name, plot_prefix, batch_key) {
  if (!"Phase" %in% colnames(obj@meta.data)) {
    warning("'Phase' column not found in metadata. Skipping plot by phase.")
    p2 <- ggplot() + theme_void() + ggtitle("UMAP by Phase (skipped)")
  } else {
    p2 <- DimPlot(obj, reduction = reduction_name, group.by = "Phase") + ggtitle(paste(plot_prefix, "UMAP by Phase"))
  }
  p1 <- DimPlot(obj, reduction = reduction_name, group.by = batch_key) + ggtitle(paste(plot_prefix, "UMAP by Batch"))
  return(list(p1 = p1, p2 = p2))
}
main <- function() {
  tryCatch({
    log_msg("Loading standardized Seurat object...")
    seurat_obj <- readRDS(input_seurat)
    
    log_msg("Performing parameter validation...")
    validation_result <- validate_params(seurat_obj, batch_key, n_features_integrate)
    seurat_obj <- validation_result$seurat_obj
    skip_correction <- validation_result$skip_correction
    
    has_scale <- FALSE
    has_scale <- tryCatch({
      nrow(GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "scale.data")) > 0
    }, error = function(e) {
      tryCatch({
        nrow(GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), slot = "scale.data")) > 0
      }, error = function(e2) FALSE)
    })
    if (!has_scale) {
      log_msg("No scale.data detected, performing ScaleData first.")
      seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    }
    log_msg("Calculating PCA and UMAP before correction...")
    seurat_pre <- RunPCA(seurat_obj, npcs = n_pcs, verbose = FALSE,seed.use=seed)
    seurat_pre <- RunUMAP(seurat_pre, dims = 1:n_pcs, reduction.name = "umap.pre", verbose = FALSE,seed.use = seed)
    plots_pre <- generate_plots(seurat_pre, "umap.pre", "Before Correction", batch_key)
    if (skip_correction) {
      log_msg("Skipping batch correction workflow.")
      seurat_corrected <- seurat_pre
      plots_post <- plots_pre
      plots_post$p1 <- plots_post$p1 + ggtitle("Correction Skipped - UMAP by Batch")
      
    } else {
      log_msg(paste("Starting batch correction workflow, method:", method))
      seurat_corrected <- run_correction(seurat_pre, method, batch_key, n_pcs, n_features_integrate)
      
      reduction_post_name <- switch(tolower(method),
                                    "harmony" = "umap.harmony",
                                    "seurat" = "umap.integrated",
                                    "combat" = "umap.combat")
      
      log_msg("Generating corrected visualization results...")
      plots_post <- generate_plots(seurat_corrected, reduction_post_name, "After Correction", batch_key)
    }
    
    log_msg("Saving visualization PDF...")
    sapply(c(output_umap, output_pca, output_corrected), function(f) {
      d <- dirname(f); if(!dir.exists(d)) dir.create(d, recursive=TRUE)
    })
    umap_layout <- (plots_pre$p1 | plots_post$p1) / (plots_pre$p2 | plots_post$p2)
    ggsave(output_umap, umap_layout, width = 12, height = 10)
    pca_plot <- DimPlot(seurat_pre, reduction = "pca", group.by = batch_key) + ggtitle("PCA Before Correction")
    ggsave(output_pca, pca_plot, width = 7, height = 6)
    log_msg("Saving corrected Seurat object...")
    saveRDS(seurat_corrected, file = output_corrected, compress = FALSE)
    
    log_msg("The process completed successfully!")
    
  }, error = function(e) {
    log_msg(sprintf("Process failed: %s", e$message), "ERROR")
    stop(e)
  })
}
main()
