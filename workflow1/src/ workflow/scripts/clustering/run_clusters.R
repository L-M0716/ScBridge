#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tools)
})


input_file       <- snakemake@input[["purified_seurat"]]
output_seurat    <- snakemake@output[["clustered_seurat"]]
output_plot      <- snakemake@output[["umap_plot_preview"]]
n_pcs            <- snakemake@params[["n_pcs"]]
reduction_method <- tolower(snakemake@params[["reduction_method"]])
n_neighbors      <- snakemake@params[["n_neighbors"]]
resolution       <- snakemake@params[["resolution"]]
algorithm_name   <- snakemake@params[["algorithm"]]
min_dist         <- snakemake@params[["min_dist"]]
spread           <- snakemake@params[["spread"]]
seed             <- snakemake@params[["seed"]]
log_file         <- snakemake@log[[1]]
sample_sheet_file <- snakemake@input[["sample_table"]]

if (!dir.exists(dirname(log_file))) dir.create(dirname(log_file), recursive = TRUE)
log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] %s", timestamp, level, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  message(line)
}

read_sample_sheet <- function(f) {
  ext <- tolower(file_ext(f))
  if (ext == "csv") {
    read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  }
}
pick_first_existing <- function(candidates, cols) {
  cand <- candidates[candidates %in% cols]
  if (length(cand) == 0) return(NA_character_) else return(cand[[1]])
}


set.seed(seed)
tryCatch({
  log_msg("Script starts execution (core clustering process).")
  
  log_msg(paste("Load Seurat object:", basename(input_file)))
  seurat_obj <- readRDS(input_file)
  
  log_msg(sprintf("Read sample sheet: %s", sample_sheet_file))
  samples_df <- read_sample_sheet(sample_sheet_file)
  sheet_id_col <- pick_first_existing(c("sample_id","Sample","sample","orig.ident","orig_ident","ID","id"),
                                      colnames(samples_df))
  group_col    <- pick_first_existing(c("group","Group","condition","Condition","treatment","Treatment"),
                                      colnames(samples_df))
  if (is.na(sheet_id_col) || is.na(group_col)) {
    stop("The sample table is missing required columns, please make sure they are included 'sample_id' and 'group'。")
  }
  log_msg(sprintf("Sample table column mapping: Sample ID column='%s'，group column='%s'。", sheet_id_col, group_col))
  
  meta_cols <- colnames(seurat_obj@meta.data)
  meta_sample_col <- pick_first_existing(c("sample_id","orig.ident","Sample","sample","orig_ident","library","Library"),
                                         meta_cols)
  if (is.na(meta_sample_col)) {
    stop("Sample ID column not found in the Seurat object's meta.data (e.g., 'sample_id' or 'orig.ident').")
  }
  log_msg(sprintf("Sample ID column in Seurat object identified as: '%s'.", meta_sample_col))
  
  seurat_obj$group <- samples_df[[group_col]][
    match(seurat_obj@meta.data[[meta_sample_col]], samples_df[[sheet_id_col]])
  ]
  group_levels <- unique(samples_df[[group_col]])
  seurat_obj$group <- factor(seurat_obj$group, levels = group_levels)
  
  uniq_meta_samples <- unique(seurat_obj@meta.data[[meta_sample_col]])
  missing_samples <- setdiff(uniq_meta_samples, samples_df[[sheet_id_col]])
  na_cells <- sum(is.na(seurat_obj$group))
  if (length(missing_samples) > 0) {
    log_msg(sprintf("WARN: The following samples are not found in the sample sheet, their group will be NA: %s",
                    paste(missing_samples, collapse = ", ")), level = "WARN")
  }
  if (na_cells > 0) {
    log_msg(sprintf("WARN: There are %d cells with NA group.", na_cells), level = "WARN")
  }
  log_msg(sprintf("Group mapping completed, detected groups: %s", paste(levels(seurat_obj$group), collapse = ", ")))
  
  reduction_for_clustering <- if (reduction_method == "harmony" && "harmony" %in% names(seurat_obj@reductions)) {
    log_msg("Harmony results detected, will use 'harmony' for clustering.")
    "harmony"
  } else {
    log_msg("No Harmony results found, will run/use 'pca' for clustering.")
    if (!"pca" %in% names(seurat_obj@reductions) || ncol(seurat_obj@reductions$pca) < n_pcs) {
      log_msg(paste("Running PCA, computing", n_pcs, "principal components..."))
      seurat_obj <- RunPCA(seurat_obj, npcs = n_pcs, verbose = FALSE)
    }
    "pca"
  }
  
  seurat_obj <- FindNeighbors(
    seurat_obj,
    reduction = reduction_for_clustering,
    dims = 1:n_pcs,
    k.param = n_neighbors,
    verbose = FALSE
  )
  
  algorithm_num <- switch(tolower(algorithm_name), "leiden" = 4, "louvain" = 1, 1)
  log_msg(sprintf("Executing %s clustering, resolution: %f", algorithm_name, resolution))
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution = resolution,
    algorithm = algorithm_num,
    verbose = FALSE
  )
  log_msg(paste("Clustering completed, found", length(unique(seurat_obj$seurat_clusters)), "clusters."))
  
  log_msg(sprintf("Running UMAP dimensionality reduction, using the first %d dimensions of '%s'.", n_pcs, reduction_for_clustering))
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = reduction_for_clustering,
    dims = 1:n_pcs,
    min.dist = min_dist,
    spread = spread,
    seed.use = seed,
    verbose = FALSE
  )
  
  log_msg("Generating cluster result preview plots (first page is the aggregate view, followed by one page per group)...")
  
  Idents(seurat_obj) <- "seurat_clusters"
  cluster_levels <- levels(Idents(seurat_obj))
  cluster_cols <- setNames(hue_pal()(length(cluster_levels)), cluster_levels)
  
  p_clusters <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, cols = cluster_cols) + 
    labs(title = "UMAP colored by Clusters") + NoLegend()
  p_batch    <- DimPlot(seurat_obj, reduction = "umap", group.by = "batch") +
    labs(title = "UMAP colored by Batch")
  page1 <- p_clusters | p_batch
  
  log_msg("Saving final clustering object and multi-page preview PDF...")
  sapply(c(output_seurat, output_plot), function(f) {
    d <- dirname(f); if(!dir.exists(d)) dir.create(d, recursive=TRUE)
  })
  saveRDS(seurat_obj, file = output_seurat, compress = FALSE)
  
  pdf(output_plot, width = 12, height = 6, onefile = TRUE)
  print(page1)
  
  groups <- levels(seurat_obj$group)
  for (g in groups) {
    if (is.na(g)) next
    log_msg(sprintf("Generating group page: %s", g))
    sub_obj <- subset(seurat_obj, subset = group == g)
    if (ncol(sub_obj) == 0) {
      log_msg(sprintf("WARN: Group '%s' contains no cells, skipping this page.", g), level = "WARN")
      next
    }
    Idents(sub_obj) <- "seurat_clusters"
    pg <- DimPlot(sub_obj, reduction = "umap", label = TRUE, repel = TRUE, cols = cluster_cols) + 
      labs(title = paste0("UMAP - Clusters (Group: ", g, ")")) + NoLegend()
    print(pg)
  }
  dev.off()
  
  log_msg("Clustering workflow completed successfully!")
  
}, error = function(e) {
  log_msg(paste("ERROR:", e$message), level = "ERROR")
  stop(e)
})