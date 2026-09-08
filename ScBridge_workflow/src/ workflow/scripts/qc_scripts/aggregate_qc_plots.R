# scripts/qc_scripts/aggregate_qc_plots.R

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(viridis)
})

seurat_qc_files <- snakemake@input[["seurat_qc_list"]]
aggregated_report_output <- snakemake@output[["aggregated_report"]]
log_file         <- snakemake@log[[1]]
if (!dir.exists(dirname(log_file))) dir.create(dirname(log_file), recursive = TRUE)
log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] [aggregate_qc] %s", timestamp, level, msg)
  cat(line, "\n", file = log_file, append = TRUE)
}

tryCatch({
  log_msg("Script started (Aggregating QC plots).")
  
  if (length(seurat_qc_files) == 0) {
    stop("No input seurat_qc.rds files found.")
  }
  
  log_msg(paste("Reading and merging", length(seurat_qc_files), "Seurat objects..."))
  
  seurat_list <- lapply(seurat_qc_files, readRDS)

  sample_ids <- gsub("/seurat_qc.rds", "", basename(dirname(seurat_qc_files)))
  names(seurat_list) <- sample_ids
  
  if (length(seurat_list) == 1) {
    merged_seurat <- seurat_list[[1]]
  } else {

    merged_seurat <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)], add.cell.ids = names(seurat_list))
  }
  
  log_msg(paste("All objects merged, total included", ncol(merged_seurat), "cells"))
  
  combined_metrics <- merged_seurat@meta.data
  log_msg("Generating aggregated QC visualization report...")
  
  p1 <- ggplot(combined_metrics, aes(x = sample_id, y = nFeature_RNA, fill = sample_id)) +
    geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
    labs(title = "Gene Count Distribution by Sample", 
         x = "Sample", y = "Number of Genes", fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position = "none")
  
  p2 <- ggplot(combined_metrics, aes(x = sample_id, y = percent.mt, fill = sample_id)) +
    geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
    labs(title = "Mitochondrial Percentage by Sample", 
         x = "Sample", y = "Mitochondrial Percentage (%)", fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position = "none")
  
  p3 <- ggplot(combined_metrics, aes(x = nCount_RNA, y = nFeature_RNA, color = sample_id)) +
    geom_point(alpha = 0.5, size = 0.8) +
    labs(title = "Genes vs UMI Counts by Sample", 
         x = "UMI Counts", y = "Gene Counts", color = "Sample") +
    theme_minimal(base_size = 12) +
    guides(color = guide_legend(override.aes = list(size=3)))
  
  p4 <- ggplot(combined_metrics, aes(x = sample_id, fill = qc_status)) +
    geom_bar(position = "fill", width = 0.8) +
    scale_fill_manual(values = c("pass" = "#1f77b4", "fail" = "#ff7f0e")) +
    labs(title = "QC Pass Rate by Sample", 
         x = "Sample", y = "Proportion", fill = "QC Status") +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  
  pdf(aggregated_report_output, width = 14, height = 10)
  print(p1 + p2)
  print(p3)
  print(p4)
  dev.off()
  Sys.sleep(1)
  system2("touch", shQuote(aggregated_report_output))
  
  log_msg("Script executed successfully.")
  
}, error = function(e) {
  log_msg(paste("ERROR:", e$message), level = "ERROR")
  stop(e)
})
