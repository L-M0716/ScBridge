# scripts/filtering_scripts/filter_cells_one_sample.R

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

input_file      <- snakemake@input[["seurat_with_doublets"]]
output_seurat   <- snakemake@output[["filtered_seurat"]]
output_stats    <- snakemake@output[["filter_stats"]]
output_report   <- snakemake@output[["filter_report"]]
min_genes       <- snakemake@params[["min_genes"]]
max_genes       <- snakemake@params[["max_genes"]]
min_counts      <- snakemake@params[["min_counts"]]
max_mito        <- snakemake@params[["max_mito"]]
min_cells_per_gene <- snakemake@params[["min_cells_per_gene"]]
exclude_genes   <- snakemake@params[["exclude_genes"]]
exclude_genes <- if (is.null(exclude_genes) || exclude_genes == "") {
  character(0)
} else {
  trimws(strsplit(exclude_genes, ",", fixed = TRUE)[[1]])
}
log_file        <- snakemake@log[[1]]
sample_id       <- snakemake@wildcards[["sample_id"]]

if (!dir.exists(dirname(log_file))) dir.create(dirname(log_file), recursive = TRUE)
log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] [%s] %s", timestamp, level, sample_id, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  message(line)
}

tryCatch({
  log_msg("Script started (cell and gene filtering).")
  
  log_msg(paste("Load Seurat object with dual-cell markers:", basename(input_file)))
  seurat_obj <- readRDS(input_file)
  initial_cells <- ncol(seurat_obj)
  log_msg(paste("Initial cell count:", initial_cells))
  if (!"qc_status" %in% colnames(seurat_obj@meta.data) || !"doublet_call" %in% colnames(seurat_obj@meta.data)) {
    stop("Input object is missing 'qc_status' or 'doublet_call' column. Please check the upstream steps.")
  }
  
  log_msg("Filtering cells based on QC status and doublet markers...")
  cells_to_keep <- rownames(seurat_obj@meta.data)[seurat_obj$qc_status == "pass" & seurat_obj$doublet_call == 0]
  
  if (length(cells_to_keep) == 0) {
    stop("No cells passed the filtering criteria. Please check the threshold settings.")
  }
  seurat_filtered <- subset(seurat_obj, cells = cells_to_keep)
  filtered_cells <- ncol(seurat_filtered)
  log_msg(sprintf("Cell filtering completed: %d / %d cells retained (%.1f%%)", 
                  filtered_cells, initial_cells, (filtered_cells/initial_cells)*100))
  initial_genes <- nrow(seurat_filtered)
  log_msg("Filtering genes...")
  
  counts_matrix <- GetAssayData(seurat_filtered, assay = "RNA", layer = "counts")
  genes_to_keep_by_expr <- Matrix::rowSums(counts_matrix > 0) >= min_cells_per_gene
  
  exclude_pattern <- paste(exclude_genes, collapse = "|")
  genes_to_keep_by_name <- !grepl(exclude_pattern, rownames(seurat_filtered), ignore.case = TRUE)
  
  final_genes_to_keep <- genes_to_keep_by_expr & genes_to_keep_by_name
  seurat_filtered <- seurat_filtered[final_genes_to_keep, ]
  filtered_genes <- nrow(seurat_filtered)
  log_msg(sprintf("Gene filtering completed: %d / %d genes retained", filtered_genes, initial_genes))
  
  stats_summary <- data.frame(
    category = c("Cells_Initial", "Cells_Kept", "Cells_Removed", 
                 "Genes_Initial (after cell filter)", "Genes_Kept", "Genes_Removed"),
    count = c(initial_cells, filtered_cells, initial_cells - filtered_cells,
              initial_genes, filtered_genes, initial_genes - filtered_genes)
  )
  
  log_msg("Generating filtering visualization report...")
  
  filter_status_df <- data.frame(
    cell_id      = colnames(seurat_obj),
    nCount_RNA   = as.numeric(seurat_obj$nCount_RNA),
    nFeature_RNA = as.numeric(seurat_obj$nFeature_RNA),
    percent_mt   = as.numeric(seurat_obj$percent.mt),
    status       = ifelse(colnames(seurat_obj) %in% colnames(seurat_filtered), "retained", "filtered"),
    stringsAsFactors = FALSE
  )
  
  p1 <- ggplot(filter_status_df, aes(x = nFeature_RNA, fill = status)) +
    geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = c(min_genes, max_genes), linetype = "dashed", color = "red") +
    labs(title = "Gene Count Distribution", x = "Number of Genes (log10 scale)", y = "Cell Count") +
    theme_bw() +
    scale_x_log10(labels = scales::label_number()) +
    annotation_logticks(sides = "b")
  
  p2 <- ggplot(filter_status_df, aes(x = percent_mt, fill = status)) +
    geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = max_mito, linetype = "dashed", color = "red") +
    labs(title = "Mitochondrial Percentage Distribution", x = "Mitochondrial %", y = "Cell Count") +
    theme_bw()
  
  p3 <- ggplot(filter_status_df, aes(x = nCount_RNA, y = nFeature_RNA, color = status)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c("retained" = "black", "filtered" = "red")) +
    geom_vline(xintercept = min_counts, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = c(min_genes, max_genes), linetype = "dashed", color = "blue") +
    labs(title = "Genes vs UMI Counts", x = "UMI Count (log10 scale)", y = "Number of Genes (log10 scale)") +
    theme_bw() +
    scale_x_log10(labels = scales::label_number()) +
    scale_y_log10(labels = scales::label_number()) +
    annotation_logticks(sides = "bl")
  
  combined_plot <- (p1 + p2) / p3 + 
    plot_annotation(title = paste("Filtering Report for Sample:", sample_id),
                    theme = theme(plot.title = element_text(hjust = 0.5, size=16)))
  
  sapply(c(output_seurat, output_stats, output_report), function(f) {
    d <- dirname(f); if(!dir.exists(d)) dir.create(d, recursive=TRUE)
  })
  
  saveRDS(seurat_filtered, file = output_seurat)
  log_msg(paste("Filtered Seurat object saved to:", basename(output_seurat)))
  
  write.csv(stats_summary, file = output_stats, row.names = FALSE)
  log_msg(paste("Filtering statistics saved to:", basename(output_stats)))
  
  ggsave(output_report, plot = combined_plot, width = 12, height = 8)
  log_msg(paste("Filtering visualization report saved to:", basename(output_report)))
  
  log_msg("Script executed successfully.")
  
}, error = function(e) {
  log_msg(paste("ERROR:", e$message), level = "ERROR")
  stop(e)
})