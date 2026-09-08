# scripts/cellchat/run_cellchat_viz_network.R

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
  library(tidyverse)
  library(igraph)
  library(ComplexHeatmap)
  library(circlize)
})

input_cellchat <- snakemake@input[["cellchat_obj"]]
output_plot_dir <- snakemake@output[["plot_dir"]]
output_circle <- snakemake@output[["circle_plot"]]
output_heatmap_count <- snakemake@output[["heatmap_count"]]
output_heatmap_weight <- snakemake@output[["heatmap_weight"]]
log_file <- snakemake@log[[1]]

if (!dir.exists(dirname(log_file))) dir.create(dirname(log_file), recursive = TRUE)
if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] %s", timestamp, level, msg)
  cat(line, "\n", file = log_file, append = TRUE)
}

write_placeholder_pdf <- function(path, text) {
  pdf(path, width = 10, height = 6)
  plot.new()
  text(0.5, 0.5, text, cex = 1.2)
  dev.off()
}

is_merged_cellchat <- function(cellchat) {
  return(!is.null(cellchat@meta$datasets) || 
         (is.list(cellchat@net) && !is.null(names(cellchat@net)) && 
          length(names(cellchat@net)) > 2))
}

tryCatch({
  
  log_msg("========== CellChat Network visualization starts ==========")
  
  if (!dir.exists(output_plot_dir)) {
    dir.create(output_plot_dir, recursive = TRUE)
  }
  
  log_msg("Loading CellChat object...")
  cellchat <- readRDS(input_cellchat)
  log_msg("CellChat object loaded successfully")
  
  if (is.list(cellchat) && isTRUE(cellchat$skipped)) {
    msg <- ifelse(is.null(cellchat$reason), "CellChat step skipped", cellchat$reason)
    if (!dir.exists(output_plot_dir)) dir.create(output_plot_dir, recursive = TRUE)
    write_placeholder_pdf(output_circle, msg)
    write_placeholder_pdf(output_heatmap_count, msg)
    write_placeholder_pdf(output_heatmap_weight, msg)
    log_msg(paste("Detection of placeholder CellChat object, placeholder plot generated:", msg), level = "WARN")
    quit(save = "no", status = 0)
  }
  
  is_merged <- is_merged_cellchat(cellchat)
  
  cellchat_list_file <- file.path(dirname(input_cellchat), "cellchat_list.rds")
  has_list <- file.exists(cellchat_list_file)
  
  if (has_list) {
    log_msg("Detection of grouped CellChat list, group comparison visualization will be performed...")
    cellchat_list <- readRDS(cellchat_list_file)
    group_names <- names(cellchat_list)
    log_msg(paste("Groups:", paste(group_names, collapse = ", ")))
    
    
    log_msg("Generating group comparison plots...")
    
    pdf(output_circle, width = 16, height = 8)
    
    par(mfrow = c(1, 2))
    
    for (i in seq_along(cellchat_list)) {
      grp <- names(cellchat_list)[i]
      cc <- cellchat_list[[grp]]
      groupSize <- as.numeric(table(cc@idents))
      
      netVisual_circle(cc@net$count, 
                       vertex.weight = groupSize, 
                       weight.scale = TRUE, 
                       label.edge = FALSE, 
                       title.name = paste(grp, "- Number of Interactions"),
                       vertex.label.cex = 0.7)
    }
    
    par(mfrow = c(1, 2))
    for (i in seq_along(cellchat_list)) {
      grp <- names(cellchat_list)[i]
      cc <- cellchat_list[[grp]]
      groupSize <- as.numeric(table(cc@idents))
      
      netVisual_circle(cc@net$weight, 
                       vertex.weight = groupSize, 
                       weight.scale = TRUE, 
                       label.edge = FALSE, 
                       title.name = paste(grp, "- Interaction Strength"),
                       vertex.label.cex = 0.7)
    }
    
    dev.off()
    log_msg("Group comparison circle plots saved")
    
    if (length(cellchat_list) == 2) {
      log_msg("Generating two-group differential comparison plots...")
      
      pdf(file.path(output_plot_dir, "comparison_differential.pdf"), 
          width = 14, height = 10)
      
      tryCatch({
        gg1 <- compareInteractions(cellchat, show.legend = TRUE, group = c(1, 2))
        gg2 <- compareInteractions(cellchat, show.legend = TRUE, group = c(1, 2), measure = "weight")
        print(gg1 + gg2)
        log_msg("Number of Interactions/Interaction Strength Comparison Plots Generated Successfully")
      }, error = function(e) {
        log_msg(paste("Interaction Comparison Plots Generation Failed:", e$message), level = "WARN")
      })
      
      tryCatch({
        par(mfrow = c(1, 2))
        netVisual_diffInteraction(cellchat, weight.scale = TRUE)
        netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
        log_msg("Differential network plots generated successfully")
      }, error = function(e) {
        log_msg(paste("Differential network plots generation failed:", e$message), level = "WARN")
      })
      
      tryCatch({
        gg1 <- netVisual_heatmap(cellchat)
        gg2 <- netVisual_heatmap(cellchat, measure = "weight")
        print(gg1 + gg2)
        log_msg("Differential heatmaps generated successfully")
      }, error = function(e) {
        log_msg(paste("Differential heatmaps generation failed:", e$message), level = "WARN")
      })
      
      dev.off()
    }
    
    log_msg("Generating heatmaps for each group...")
    
    pdf(output_heatmap_count, width = 12, height = 10)
    for (grp in names(cellchat_list)) {
      tryCatch({
        netVisual_heatmap(cellchat_list[[grp]], 
                          measure = "count",
                          title.name = paste(grp, "- Number of Interactions"))
      }, error = function(e) {
        log_msg(paste("Group", grp, "heatmap generation failed:", e$message), level = "WARN")
      })
    }
    dev.off()
    
    pdf(output_heatmap_weight, width = 12, height = 10)
    for (grp in names(cellchat_list)) {
      tryCatch({
        netVisual_heatmap(cellchat_list[[grp]], 
                          measure = "weight",
                          title.name = paste(grp, "- Interaction Strength"))
      }, error = function(e) {
        log_msg(paste("Group", grp, "heatmap generation failed:", e$message), level = "WARN")
      })
    }
    dev.off()
    
    log_msg("Group heatmaps saved")
    
  } else {
    log_msg("Performing single-group visualization...")
    
    groupSize <- as.numeric(table(cellchat@idents))
    n_celltypes <- length(groupSize)
    log_msg(paste("Detected", n_celltypes, "cell types"))
    
    pdf(output_circle, width = 14, height = 7)
    par(mfrow = c(1, 2), xpd = TRUE)
    
    netVisual_circle(cellchat@net$count, 
                     vertex.weight = groupSize, 
                     weight.scale = TRUE, 
                     label.edge = FALSE, 
                     title.name = "Number of Interactions",
                     vertex.label.cex = 0.7)
    
    netVisual_circle(cellchat@net$weight, 
                     vertex.weight = groupSize, 
                     weight.scale = TRUE, 
                     label.edge = FALSE, 
                     title.name = "Interaction Strength",
                     vertex.label.cex = 0.7)
    dev.off()
    
    pdf(output_heatmap_count, width = 10, height = 8)
    netVisual_heatmap(cellchat, 
                      measure = "count",
                      color.heatmap = c("#2166ac", "#f7f7f7", "#b2182b"),
                      title.name = "Number of Interactions")
    dev.off()
    
    pdf(output_heatmap_weight, width = 10, height = 8)
    netVisual_heatmap(cellchat, 
                      measure = "weight",
                      color.heatmap = c("#2166ac", "#f7f7f7", "#b2182b"),
                      title.name = "Interaction Strength")
    dev.off()
  }
  
  
  log_msg("Generating bubble plots...")
  
  if (has_list) {
    pdf(file.path(output_plot_dir, "interaction_bubble_comparison.pdf"), 
        width = 16, height = 12)
    
    for (grp in names(cellchat_list)) {
      cc <- cellchat_list[[grp]]
      n_celltypes <- length(levels(cc@idents))
      
      tryCatch({
        p <- netVisual_bubble(cc, 
                              sources.use = 1:min(10, n_celltypes), 
                              targets.use = 1:n_celltypes, 
                              remove.isolate = TRUE) +
          ggtitle(paste("Group:", grp))
        print(p)
      }, error = function(e) {
        log_msg(paste("Group", grp, "bubble plot generation failed:", e$message), level = "WARN")
      })
    }
    dev.off()
  } else {
    pdf(file.path(output_plot_dir, "interaction_bubble_all.pdf"), 
        width = 14, height = max(10, n_celltypes * 0.5))
    
    tryCatch({
      p <- netVisual_bubble(cellchat, 
                            sources.use = 1:min(15, n_celltypes), 
                            targets.use = 1:n_celltypes, 
                            remove.isolate = TRUE)
      print(p)
    }, error = function(e) {
      log_msg(paste("Bubble plot generation failed:", e$message), level = "WARN")
    })
    dev.off()
  }
  
  log_msg("========== CellChat Network Visualization Complete ==========")
  
}, error = function(e) {
  log_msg(paste("Error:", e$message), level = "ERROR")
  stop(e)
})
