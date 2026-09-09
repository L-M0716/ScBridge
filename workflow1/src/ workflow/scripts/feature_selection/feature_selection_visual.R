# scripts/feature_selection_scripts/visualize_and_purify_features.R
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

input_seurat   <- snakemake@input[["corrected_seurat"]]
output_seurat <- snakemake@output[["feature_seurat"]]
output_plot    <- snakemake@output[["feature_plot"]]
log_file       <- snakemake@log[[1]]

if (!dir.exists(dirname(log_file))) dir.create(dirname(log_file), recursive = TRUE)
log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] %s", timestamp, level, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  message(line)
}

tryCatch({
  log_msg("Script execution started (highly variable gene visualization).")
  
  log_msg(paste("Load the batch-corrected Seurat object:", basename(input_seurat)))
  seurat_obj <- readRDS(input_seurat)
  
  available_assays <- Assays(seurat_obj)
  log_msg(paste("Available assays:", paste(available_assays, collapse = ", ")))
  
  original_assay <- DefaultAssay(seurat_obj)
  
  if ("integrated" %in% available_assays) {
    DefaultAssay(seurat_obj) <- "integrated"
    hvg_integrated <- VariableFeatures(seurat_obj)
    log_msg(paste("从integrated assay获取", length(hvg_integrated), "个高变基因"))
  } else {
    hvg_integrated <- NULL
  }
  
  viz_assay <- if ("SCT" %in% available_assays) {
    "SCT"
  } else if ("RNA" %in% available_assays) {
    "RNA"
  } else {
    original_assay
  }
  
  DefaultAssay(seurat_obj) <- viz_assay
  log_msg(paste("use", viz_assay, "Visualizing the assay"))
  

  if (!is.null(hvg_integrated) && length(hvg_integrated) > 0) {

    hvg <- hvg_integrated
    VariableFeatures(seurat_obj) <- hvg
    log_msg("List of highly variable genes using the integrated assay")
  } else {
    hvg <- VariableFeatures(seurat_obj)
    if (length(hvg) == 0) {
      log_msg("No highly variable genes found, recalculating...")
      seurat_obj <- FindVariableFeatures(seurat_obj, 
                                         selection.method = "vst", 
                                         nfeatures = 2000,
                                         verbose = FALSE)
      hvg <- VariableFeatures(seurat_obj)
    }
  }
  
  log_msg(paste("Number of highly variable genes for visualization:", length(hvg)))

  log_msg("Recalculating HVF information for visualization...")
  seurat_obj <- FindVariableFeatures(seurat_obj, 
                                     selection.method = "vst", 
                                     nfeatures = min(3000, nrow(seurat_obj)),
                                     verbose = FALSE)
  
  log_msg("Preparing visual charts...")
  
  p1 <- tryCatch({
    VariableFeaturePlot(seurat_obj) +
      labs(title = paste("Highly Variable Genes (", viz_assay, " assay)", sep = "")) +
      theme_minimal()
  }, error = function(e) {
    log_msg(paste("VariableFeaturePlot failed:", e$message), "WARN")
    NULL
  })
  
  p2 <- tryCatch({
    hvf_info <- HVFInfo(seurat_obj, selection.method = "vst")
    
    if (is.null(hvf_info) || nrow(hvf_info) == 0) {
      stop("HVF information is empty")
    }
    
    hvf_info$gene <- rownames(hvf_info)
    
    col_names <- colnames(hvf_info)
    mean_col <- col_names[grep("mean|avg", col_names, ignore.case = TRUE)][1]
    var_col <- col_names[grep("variance|dispersion|residual_variance", col_names, ignore.case = TRUE)][1]
    
    if (is.na(mean_col)) mean_col <- "mean"
    if (is.na(var_col)) var_col <- "variance.standardized"
    
    log_msg(paste("Using columns:", mean_col, "and", var_col))
    

    hvf_info$is_hvg <- hvf_info$gene %in% hvg
    

    top10_hvg <- hvg[hvg %in% hvf_info$gene][1:min(10, sum(hvg %in% hvf_info$gene))]
    

    p <- ggplot(hvf_info, aes_string(x = mean_col, y = var_col)) +
      geom_point(aes(color = is_hvg), alpha = 0.5, size = 0.8) +
      scale_color_manual(values = c("TRUE" = "tomato", "FALSE" = "grey80"),
                         labels = c("TRUE" = "HVG", "FALSE" = "Not HVG")) +
      labs(
        title = "Top Highly Variable Genes",
        subtitle = paste(sum(hvf_info$is_hvg), "genes identified as variable"),
        x = "Average Expression", 
        y = "Standardized Variance",
        color = "Gene Type"
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom")
    
   
    if (length(top10_hvg) > 0) {
      top10_data <- hvf_info[hvf_info$gene %in% top10_hvg, ]
      p <- p + geom_text_repel(
        data = top10_data,
        aes(label = gene), 
        color = "black", 
        size = 3.5,
        box.padding = 0.5, 
        max.overlaps = 20
      )
    }
    
    p
    
  }, error = function(e) {
    log_msg(paste("Failed to generate HVG scatter plot:", e$message), "WARN")
    NULL
  })
  

  p3 <- tryCatch({

    top_genes <- head(hvg, 30)
    

    expr_data <- GetAssayData(seurat_obj, slot = "data")
    mean_expr <- rowMeans(expr_data[top_genes, , drop = FALSE])
    
    df <- data.frame(
      gene = names(mean_expr),
      mean_expression = mean_expr
    )
    df <- df[order(df$mean_expression, decreasing = TRUE), ]
    df$gene <- factor(df$gene, levels = df$gene)
    
    ggplot(df[1:min(20, nrow(df)), ], aes(x = gene, y = mean_expression)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(
        title = "Top 20 Highly Variable Genes",
        x = "Gene",
        y = "Mean Expression"
      ) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8))
    
  }, error = function(e) {
    log_msg(paste("Failed to generate bar chart:", e$message), "WARN")
    NULL
  })
  

  p4 <- NULL
  available_reductions <- names(seurat_obj@reductions)
  log_msg(paste("Available dimensionality reduction results:", paste(available_reductions, collapse = ", ")))
  
  reduction_to_use <- if ("umap.integrated" %in% available_reductions) {
    "umap.integrated"
  } else if ("umap.harmony" %in% available_reductions) {
    "umap.harmony"
  } else if ("umap" %in% available_reductions) {
    "umap"
  } else {
    NULL
  }
  
  if (!is.null(reduction_to_use)) {
    log_msg(paste("Using dimensionality reduction results:", reduction_to_use))
    

    if ("integrated" %in% available_assays) {
      DefaultAssay(seurat_obj) <- "integrated"
    }
    

    available_features <- rownames(seurat_obj)
    display_genes <- hvg[hvg %in% available_features][1:min(12, sum(hvg %in% available_features))]
    
    if (length(display_genes) > 0) {
      p4 <- FeaturePlot(seurat_obj, 
                        features = display_genes,
                        reduction = reduction_to_use,
                        ncol = 4, 
                        pt.size = 0.1) & 
        theme(axis.title = element_blank(), 
              axis.text = element_blank(),
              axis.ticks = element_blank(), 
              plot.title = element_text(size = 10))
    }
  }
  
  log_msg("Combining charts and saving PDF report...")

  sapply(c(output_seurat, output_plot), function(f) {
    d <- dirname(f)
    if(!dir.exists(d)) dir.create(d, recursive = TRUE)
  })
  
  pdf(output_plot, width = 12, height = 14)
  
  plots_available <- list(p3, p4)
  plots_valid <- plots_available[!sapply(plots_available, is.null)]
  
  if (length(plots_valid) > 0) {
    if (length(plots_valid) == 1) {
      print(plots_valid[[1]])
    } else {
      print(plots_valid[[1]] / plots_valid[[2]])
    }
  } else {
    print(ggplot() + 
            annotate("text", x = 0.5, y = 0.5, 
                     label = "No visualization generated due to errors", 
                     size = 6) + 
            theme_void())
  }
  
  dev.off()
  log_msg(paste("Visualization report successfully saved to:", basename(output_plot)))
  

  DefaultAssay(seurat_obj) <- original_assay
  

  log_msg(paste("Saving Seurat object to:", basename(output_seurat)))
  saveRDS(seurat_obj, file = output_seurat, compress = FALSE)
  

  log_msg("========== Highly Variable Genes Statistics ==========")
  log_msg(paste("Total number of highly variable genes:", length(hvg)))
  log_msg(paste("Top 10 highly variable genes:", paste(head(hvg, 10), collapse = ", ")))
  log_msg(paste("Current default assay:", DefaultAssay(seurat_obj)))
  log_msg("===================================")
  
  log_msg("Script executed successfully.")
  
}, error = function(e) {
  log_msg(paste("ERROR:", e$message), level = "ERROR")
  

  tryCatch({
    log_msg("Trying to save the original input object...")
    seurat_obj <- readRDS(input_seurat)
    saveRDS(seurat_obj, file = output_seurat, compress = FALSE)
    log_msg("Original object saved.")
  }, error = function(e2) {
    log_msg(paste("Failed to save object:", e2$message), "ERROR")
  })
  
  stop(e)
})