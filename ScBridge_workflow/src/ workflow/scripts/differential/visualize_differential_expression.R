# scripts/differential_scripts/visualize_differential_expression.R

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(patchwork)
})
de_results_file   <- snakemake@input[["de_results_table"]]
seurat_obj_path   <- snakemake@input[["annotated_seurat"]]
output_volcano    <- snakemake@output[["volcano_plots"]]
output_heatmap    <- snakemake@output[["de_heatmap"]]
output_summary    <- snakemake@output[["top_de_genes_table"]]
log_file          <- snakemake@log[[1]]
top_n             <- snakemake@params[["top_n"]]
adj_pval_cutoff   <- snakemake@params[["adj_pval_cutoff"]]
vis_logfc_thresh  <- snakemake@params[["vis_logfc_threshold"]]
sapply(c(output_volcano, output_heatmap, output_summary, log_file), function(f) {
  d <- dirname(f); if(!dir.exists(d)) dir.create(d, recursive=TRUE)
})
log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] %s", timestamp, level, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  message(line)
}

tryCatch({
  log_msg("Start the process of visualizing the differential analysis results.")
  
  log_msg(paste("Load difference analysis results:", basename(de_results_file)))
  de_results <- read.csv(de_results_file)
  
  if (nrow(de_results) == 0) {
    log_msg("Differential analysis results are empty, cannot generate visualization plots. Creating empty files and exiting.", level="WARN")
    file.create(output_volcano)
    file.create(output_heatmap)
    file.create(output_summary)
    quit(save = "no", status = 0)
  }
  
  log_msg("Generating volcano plots...")
  pdf(output_volcano, width = 8, height = 7)
  
  cell_types_in_res <- unique(de_results$cell_type)
  for (ct in cell_types_in_res) {
    df_subset <- de_results %>% 
      filter(cell_type == ct) %>%
      mutate(significance = ifelse(p_val_adj < adj_pval_cutoff & abs(avg_log2FC) > vis_logfc_thresh, 
                                   paste0("p.adj < ", adj_pval_cutoff, " & |logFC| > ", vis_logfc_thresh), 
                                   "Not Significant"))
    
    df_subset$label <- ifelse(df_subset$significance != "Not Significant", df_subset$gene, "")
    
    p <- ggplot(df_subset, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
      geom_point(aes(color = significance), alpha = 0.6) +
      geom_text_repel(aes(label = label), max.overlaps = 15, size = 3,
                      box.padding = 0.5, point.padding = 0.2,
                      min.segment.length = 0, seed = 42) +
      scale_color_manual(values = c("Not Significant" = "grey80",
                                   setNames("red", paste0("p.adj < ", adj_pval_cutoff,
                                                          " & |logFC| > ", vis_logfc_thresh))),
                       name = "Significance") +
      geom_vline(xintercept = c(-vis_logfc_thresh, vis_logfc_thresh), linetype = "dashed") +
      geom_hline(yintercept = -log10(adj_pval_cutoff), linetype = "dashed") +
      labs(title = paste("Volcano Plot for", ct),
           subtitle = unique(df_subset$comparison),
           x = "log2 Fold Change", y = "-log10(Adjusted P-value)") +
      theme_bw() +
      theme(legend.position = "bottom")
    
    print(p)
  }
  dev.off()
  log_msg(paste("Volcano plots saved to:", basename(output_volcano)))
  
  log_msg("Filtering top N genes for heatmap and summary table...")
  
  top_genes <- de_results %>%
    filter(p_val_adj < adj_pval_cutoff) %>%
    group_by(cell_type) %>%
    slice_max(order_by = avg_log2FC, n = ceiling(top_n / 2)) %>%
    bind_rows(
      de_results %>%
        filter(p_val_adj < adj_pval_cutoff) %>%
        group_by(cell_type) %>%
        slice_min(order_by = avg_log2FC, n = floor(top_n / 2))
    ) %>%
    ungroup() %>%
    distinct(gene, .keep_all = TRUE)
  write.csv(top_genes, file = output_summary, row.names = FALSE)
  log_msg(paste("Top differentially expressed genes summary table saved to:", basename(output_summary)))
  if (nrow(top_genes) > 0) {
    log_msg("Loading Seurat object to generate heatmap...")
    seurat_obj <- readRDS(seurat_obj_path)
    
    group_by_col <- "group"
    celltype_col <- if ("celltype_fine" %in% colnames(seurat_obj@meta.data)) {
      "celltype_fine"
    } else if ("singler_labels" %in% colnames(seurat_obj@meta.data)) {
      "singler_labels"
    } else if ("manual_annotation" %in% colnames(seurat_obj@meta.data)) {
      "manual_annotation"
    } else {
      log_msg("Warning: Cell type annotation column not found, using seurat_clusters", level = "WARN")
      "seurat_clusters"
    }
    
    if (!group_by_col %in% colnames(seurat_obj@meta.data)) {
      log_msg("Warning: 'group' column not found in Seurat object, attempting to use other grouping information", level = "WARN")
      if ("sample_id" %in% colnames(seurat_obj@meta.data)) {
        seurat_obj$celltype_group <- paste(seurat_obj@meta.data[[celltype_col]], 
                                           seurat_obj@meta.data[["sample_id"]], sep = "_")
      } else {
        seurat_obj$celltype_group <- seurat_obj@meta.data[[celltype_col]]
      }
    } else {
      seurat_obj$celltype_group <- paste(seurat_obj@meta.data[[celltype_col]], 
                                         seurat_obj@meta.data[[group_by_col]], sep = "_")
    }
    
    seurat_obj$celltype_group <- factor(seurat_obj$celltype_group, 
                                        levels = sort(unique(seurat_obj$celltype_group)))
    
    Idents(seurat_obj) <- "celltype_group"

    available_assays <- names(seurat_obj@assays)
    if ("assay_used" %in% colnames(de_results)) {
      de_assay  <- unique(de_results$assay_used)[1]
      use_assay <- if (de_assay %in% available_assays) de_assay else
                   if ("RNA" %in% available_assays) "RNA" else "SCT"
    } else {
      use_assay <- if ("RNA" %in% available_assays) "RNA" else
                   if ("SCT" %in% available_assays) "SCT" else available_assays[1]
    }
    log_msg(paste("Heatmap uses assay:", use_assay))

    DefaultAssay(seurat_obj) <- use_assay

    available_genes <- rownames(seurat_obj[[use_assay]])
    valid_genes <- intersect(top_genes$gene, available_genes)
    
    if (length(valid_genes) == 0) {
      log_msg("Warning: Top genes do not exist in the current assay, cannot generate heatmap", level = "WARN")
      file.create(output_heatmap)
    } else {
      if (length(valid_genes) < nrow(top_genes)) {
        log_msg(paste("Warning: ", nrow(top_genes) - length(valid_genes), 
                      " genes do not exist in the current assay"), level = "WARN")
      }
      genes_to_scale <- valid_genes[!valid_genes %in%
        rownames(GetAssayData(seurat_obj, assay = use_assay, layer = "scale.data"))]
      if (length(genes_to_scale) > 0) {
        log_msg(paste("Scaling", length(genes_to_scale), " DE genes..."))
        seurat_obj <- ScaleData(seurat_obj, features = genes_to_scale,
                                assay = use_assay, verbose = FALSE)
      }

      p_heatmap <- DoHeatmap(
        seurat_obj, 
        features = valid_genes,
        group.by = "celltype_group", 
        label = FALSE
      ) + labs(title = "Top Differentially Expressed Genes")
      
      ggsave(output_heatmap, plot = p_heatmap, 
             width = 14, 
             height = max(10, length(valid_genes) * 0.2), 
             limitsize = FALSE)
      log_msg(paste("Differentially expressed genes heatmap saved to:", basename(output_heatmap)))
    }
  } else {
    log_msg("No significant genes found to generate heatmap.", level = "WARN")
    file.create(output_heatmap)
  }
  log_msg("Differential expression analysis visualization completed successfully!")
}, error = function(e) {
  log_msg(paste("ERROR:", e$message), level = "ERROR")
  stop(e)
})