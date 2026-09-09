# scripts/cellchat/run_cellchat_viz_pathway.R

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
  library(tidyverse)
  library(ggplot2)
})

input_cellchat <- snakemake@input[["cellchat_obj"]]

get_output_path <- function(output_obj, name, index = NULL) {
  val <- tryCatch(output_obj[[name]], error = function(e) NULL)
  if (is.null(val) || length(val) == 0) {
    if (!is.null(index)) {
      val <- tryCatch(output_obj[[index]], error = function(e) NULL)
    }
  }
  if (is.null(val) || length(val) == 0) {
    all_names <- names(output_obj)
    if (name %in% all_names) {
      idx <- which(all_names == name)
      val <- output_obj[[idx]]
    }
  }
  return(val)
}

output_pathway_summary <- get_output_path(snakemake@output, "pathway_summary", 2)
output_signaling_role <- get_output_path(snakemake@output, "signaling_role", 3)

if (!is.null(output_signaling_role) && length(output_signaling_role) > 0) {
  output_plot_dir <- dirname(output_signaling_role)
} else if (!is.null(output_pathway_summary) && length(output_pathway_summary) > 0) {
  output_plot_dir <- file.path(dirname(dirname(output_pathway_summary)), "plots", "pathways")
} else {
  output_plot_dir <- get_output_path(snakemake@output, "plot_dir", 1)
}

pathways_of_interest <- snakemake@params[["pathways_of_interest"]]
top_n_pathways <- snakemake@params[["top_n_pathways"]]
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

analyze_pathways_single <- function(cellchat, group_name, output_dir, 
                                     pathways_of_interest, top_n_pathways) {
  log_msg(paste("=== Analysis group:", group_name, "==="))
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  pathways <- cellchat@netP$pathways
  log_msg(paste("  Detected", length(pathways), "significant signaling pathways"))
  
  if (length(pathways) == 0) {
    log_msg(paste("  Warning: Group", group_name, "did not detect significant signaling pathways"), level = "WARN")
    return(list(
      pathway_summary = data.frame(message = "No significant signaling pathways detected"),
      n_pathways = 0
    ))
  }
  
  log_msg("  Computing network centrality scores...")
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  if (!is.null(pathways_of_interest) && length(pathways_of_interest) > 0 && 
      pathways_of_interest[1] != "") {
    pathways_to_analyze <- intersect(pathways_of_interest, pathways)
    if (length(pathways_to_analyze) == 0) {
      pathways_to_analyze <- head(pathways, top_n_pathways)
    }
  } else {
    pathways_to_analyze <- head(pathways, top_n_pathways)
  }
  
  log_msg(paste("  Analyzing pathways:", paste(pathways_to_analyze, collapse = ", ")))
  
  pathway_stats <- lapply(pathways, function(p) {
    tryCatch({
      df <- subsetCommunication(cellchat, signaling = p)
      data.frame(
        group = group_name,
        pathway = p,
        n_interactions = nrow(df),
        n_source_celltypes = length(unique(df$source)),
        n_target_celltypes = length(unique(df$target)),
        total_prob = round(sum(df$prob, na.rm = TRUE), 4),
        mean_prob = round(mean(df$prob, na.rm = TRUE), 4)
      )
    }, error = function(e) {
      data.frame(
        group = group_name,
        pathway = p,
        n_interactions = NA,
        n_source_celltypes = NA,
        n_target_celltypes = NA,
        total_prob = NA,
        mean_prob = NA
      )
    })
  })
  
  pathway_summary <- do.call(rbind, pathway_stats)
  pathway_summary <- pathway_summary[order(-pathway_summary$n_interactions), ]
  
  pdf(file.path(output_dir, paste0("signaling_role_scatter_", group_name, ".pdf")), 
      width = 12, height = 10)
  
  tryCatch({
    gg1 <- netAnalysis_signalingRole_scatter(cellchat) +
      ggtitle(paste("Signaling Role Analysis -", group_name)) +
      theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))
    print(gg1)
    log_msg("  Signaling role scatter plot generated successfully")
  }, error = function(e) {
    log_msg(paste("  Failed to generate signaling role scatter plot:", e$message), level = "WARN")
    plot.new()
    text(0.5, 0.5, "Failed to generate signaling role scatter plot", cex = 1.5)
  })
  
  dev.off()
  
  pathway_detail_dir <- file.path(output_dir, "pathway_details")
  if (!dir.exists(pathway_detail_dir)) dir.create(pathway_detail_dir, recursive = TRUE)
  
  n_celltypes <- length(levels(cellchat@idents))
  
  for (pathway in pathways_to_analyze) {
    log_msg(paste("    Analyzing pathway:", pathway))
    
    tryCatch({
      pdf(file.path(pathway_detail_dir, paste0(pathway, "_", group_name, ".pdf")), 
          width = 12, height = 10)
      
      tryCatch({
        par(mfrow = c(1, 1), mar = c(1, 1, 3, 1))
        netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
        title(main = paste0(pathway, " - ", group_name), cex.main = 1.2)
      }, error = function(e) {})
      
      tryCatch({
        par(mfrow = c(1, 1), mar = c(1, 1, 3, 1))
        netVisual_aggregate(cellchat, signaling = pathway, layout = "chord")
        title(main = paste0(pathway, " Chord - ", group_name), cex.main = 1.2)
      }, error = function(e) {})
      
      tryCatch({
        par(mfrow = c(1, 1))
        ht <- netVisual_heatmap(cellchat, signaling = pathway, color.heatmap = "Reds")
        print(ht)
      }, error = function(e) {})
      
      tryCatch({
        gg_contrib <- netAnalysis_contribution(cellchat, signaling = pathway)
        print(gg_contrib + ggtitle(paste0(pathway, " L-R Contribution - ", group_name)))
      }, error = function(e) {})
      
      dev.off()
      
    }, error = function(e) {
      log_msg(paste("    Pathway", pathway, "visualization failed:", e$message), level = "WARN")
      tryCatch(dev.off(), error = function(x) {})
    })
  }
  
  pdf(file.path(output_dir, paste0("pathway_bubble_", group_name, ".pdf")), 
      width = 14, height = 10)
  
  batch_size <- 5
  n_batches <- ceiling(length(pathways_to_analyze) / batch_size)
  
  for (b in 1:n_batches) {
    start_idx <- (b - 1) * batch_size + 1
    end_idx <- min(b * batch_size, length(pathways_to_analyze))
    batch_pathways <- pathways_to_analyze[start_idx:end_idx]
    
    tryCatch({
      p <- netVisual_bubble(cellchat, signaling = batch_pathways, remove.isolate = TRUE)
      print(p + ggtitle(paste(group_name, "-", paste(batch_pathways, collapse = ", "))))
    }, error = function(e) {
      log_msg(paste("  Bubble plot batch", b, "failed:", e$message), level = "WARN")
    })
  }
  
  dev.off()
  
  pdf(file.path(output_dir, paste0("gene_expression_", group_name, ".pdf")), 
      width = 12, height = 8)
  
  for (pathway in head(pathways_to_analyze, 5)) {
    tryCatch({
      p <- plotGeneExpression(cellchat, signaling = pathway)
      print(p + ggtitle(paste0(pathway, " Gene Expression - ", group_name)))
    }, error = function(e) {})
  }
  
  dev.off()
  
  return(list(
    cellchat = cellchat,
    pathway_summary = pathway_summary,
    n_pathways = length(pathways),
    pathways_analyzed = pathways_to_analyze
  ))
}

tryCatch({
  
  log_msg("========== CellChat Signaling Pathway Analysis Started ==========")
  
  if (is.null(output_plot_dir) || length(output_plot_dir) == 0 || output_plot_dir == "") {
    stop("Failed to retrieve output_plot_dir path")
  }
  if (is.null(output_pathway_summary) || length(output_pathway_summary) == 0) {
    stop("Failed to retrieve output_pathway_summary path")
  }
  if (is.null(output_signaling_role) || length(output_signaling_role) == 0) {
    stop("Failed to retrieve output_signaling_role path")
  }
  
  if (!dir.exists(output_plot_dir)) dir.create(output_plot_dir, recursive = TRUE)
  if (!dir.exists(dirname(output_pathway_summary))) {
    dir.create(dirname(output_pathway_summary), recursive = TRUE)
  }
  
  cellchat_check <- readRDS(input_cellchat)
  if (is.list(cellchat_check) && isTRUE(cellchat_check$skipped)) {
    msg <- ifelse(is.null(cellchat_check$reason), "CellChat step skipped", cellchat_check$reason)
    write.csv(data.frame(message = msg), output_pathway_summary, row.names = FALSE)
    write_placeholder_pdf(output_signaling_role, msg)
    log_msg(paste("Detected placeholder CellChat object, placeholder pathway results output:", msg), level = "WARN")
    quit(save = "no", status = 0)
  }
  
  cellchat_list_file <- file.path(dirname(input_cellchat), "cellchat_list.rds")
  has_groups <- file.exists(cellchat_list_file)
  
  all_pathway_summary <- data.frame()
  
  if (has_groups) {
    log_msg("Detected grouped CellChat list, performing grouped pathway analysis...")
    
    cellchat_list <- readRDS(cellchat_list_file)
    group_names <- names(cellchat_list)
    log_msg(paste("group:", paste(group_names, collapse = ", ")))
    
    results_list <- list()
    
    for (grp in group_names) {
      group_output_dir <- file.path(output_plot_dir, grp)
      
      result <- analyze_pathways_single(
        cellchat = cellchat_list[[grp]],
        group_name = grp,
        output_dir = group_output_dir,
        pathways_of_interest = pathways_of_interest,
        top_n_pathways = top_n_pathways
      )
      
      results_list[[grp]] <- result
      all_pathway_summary <- rbind(all_pathway_summary, result$pathway_summary)
    }
    
    if (length(cellchat_list) >= 2) {
      log_msg("\nPerforming grouped comparison analysis...")
      
      comparison_dir <- file.path(output_plot_dir, "comparison")
      if (!dir.exists(comparison_dir)) dir.create(comparison_dir, recursive = TRUE)
      
      cellchat_merged <- readRDS(input_cellchat)
      
      pdf(file.path(comparison_dir, "pathway_ranking_comparison.pdf"), 
          width = 14, height = 10)
      
      tryCatch({
        gg1 <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
        gg2 <- rankNet(cellchat_merged, mode = "comparison", stacked = FALSE, do.stat = TRUE)
        print(gg1 + gg2)
        log_msg("Pathway ranking comparison plot generated successfully")
      }, error = function(e) {
        log_msg(paste("Pathway ranking comparison plot failed:", e$message), level = "WARN")
        plot.new()
        text(0.5, 0.5, "Failed to generate pathway ranking comparison plot", cex = 1.5)
      })
      
      dev.off()
      
      pdf(file.path(comparison_dir, "information_flow_comparison.pdf"), 
          width = 12, height = 8)
      
      tryCatch({
        gg <- rankNet(cellchat_merged, mode = "comparison", stacked = FALSE, 
                      do.stat = TRUE, return.data = FALSE)
        print(gg + ggtitle("Information Flow Comparison"))
        log_msg("Information flow comparison plot generated successfully")
      }, error = function(e) {
        log_msg(paste("Information flow comparison plot failed:", e$message), level = "WARN")
      })
      
      dev.off()
      
      if (length(group_names) == 2) {
        log_msg("Identifying significantly different pathways between groups...")
        
        pathways_1 <- results_list[[group_names[1]]]$cellchat@netP$pathways
        pathways_2 <- results_list[[group_names[2]]]$cellchat@netP$pathways
        common_pathways <- intersect(pathways_1, pathways_2)
        unique_to_1 <- setdiff(pathways_1, pathways_2)
        unique_to_2 <- setdiff(pathways_2, pathways_1)
        
        log_msg(paste("  Common pathways:", length(common_pathways)))
        log_msg(paste("  Only", group_names[1], "has:", length(unique_to_1)))
        log_msg(paste("  Only", group_names[2], "has:", length(unique_to_2)))
        
        pathway_comparison <- data.frame(
          pathway = c(common_pathways, unique_to_1, unique_to_2),
          status = c(rep("common", length(common_pathways)),
                     rep(paste("unique_to", group_names[1]), length(unique_to_1)),
                     rep(paste("unique_to", group_names[2]), length(unique_to_2)))
        )
        
        write.csv(pathway_comparison, 
                  file.path(comparison_dir, "pathway_comparison.csv"),
                  row.names = FALSE)
        
        if (length(common_pathways) > 0) {
          pdf(file.path(comparison_dir, "signaling_role_comparison.pdf"), 
              width = 16, height = 10)
          
          for (grp in group_names) {
            tryCatch({
              cc <- results_list[[grp]]$cellchat
              ht <- netAnalysis_signalingRole_heatmap(cc, pattern = "all", 
                                                       signaling = head(common_pathways, 20),
                                                       title = paste("Signaling Role -", grp))
              print(ht)
            }, error = function(e) {
              log_msg(paste("Group", grp, "signaling role heatmap failed:", e$message), level = "WARN")
            })
          }
          
          dev.off()
        }
        
        pdf(file.path(comparison_dir, "top_pathway_comparison.pdf"), 
            width = 16, height = 12)
        
        for (pathway in head(common_pathways, 5)) {
          log_msg(paste("  Comparing pathway:", pathway))
          
          tryCatch({
            par(mfrow = c(1, 2), mar = c(2, 2, 3, 2))
            
            for (grp in group_names) {
              cc <- results_list[[grp]]$cellchat
              netVisual_aggregate(cc, signaling = pathway, layout = "circle")
              title(main = paste(pathway, "-", grp), cex.main = 1.2)
            }
            
          }, error = function(e) {
            log_msg(paste("  Pathway", pathway, "comparison plot failed:", e$message), level = "WARN")
          })
        }
        
        dev.off()
      }
    }
    
    pdf(output_signaling_role, width = 14, height = 10)
    
    tryCatch({
      cellchat_merged <- readRDS(input_cellchat)
      gg <- netAnalysis_signalingRole_scatter(cellchat_merged) +
        ggtitle("Signaling Role Analysis - All Groups") +
        theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))
      print(gg)
    }, error = function(e) {
      if (length(results_list) > 0) {
        cc <- results_list[[1]]$cellchat
        gg <- netAnalysis_signalingRole_scatter(cc) +
          ggtitle(paste("Signaling Role Analysis -", names(results_list)[1]))
        print(gg)
      } else {
        plot.new()
        text(0.5, 0.5, "Failed to generate signaling role scatter plot", cex = 1.5)
      }
    })
    
    dev.off()
    
  } else {
    log_msg("Perform single-group pathway analysis...")
    
    cellchat <- readRDS(input_cellchat)
    
    result <- analyze_pathways_single(
      cellchat = cellchat,
      group_name = "all",
      output_dir = output_plot_dir,
      pathways_of_interest = pathways_of_interest,
      top_n_pathways = top_n_pathways
    )
    
    all_pathway_summary <- result$pathway_summary
    
    file.copy(
      file.path(output_plot_dir, "signaling_role_scatter_all.pdf"),
      output_signaling_role,
      overwrite = TRUE
    )
  }
  
  write.csv(all_pathway_summary, output_pathway_summary, row.names = FALSE, quote = FALSE)
  log_msg(paste("Pathway summary table saved to:", basename(output_pathway_summary)))
  
  log_msg("========== CellChat signaling pathway analysis completed ==========")
  
}, error = function(e) {
  log_msg(paste("Error:", e$message), level = "ERROR")
  log_msg(paste("Error traceback:", paste(capture.output(traceback()), collapse = "\n")), level = "ERROR")
  stop(e)
})
