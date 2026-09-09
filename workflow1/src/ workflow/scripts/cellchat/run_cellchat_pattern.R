# scripts/cellchat/run_cellchat_pattern.R
# 支持按group分组分析 + 分组比较

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
  library(tidyverse)
  library(ggplot2)
  library(NMF)
  library(ggalluvial)
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

output_cellchat_final <- get_output_path(snakemake@output, "cellchat_final", 1)
output_outgoing <- get_output_path(snakemake@output, "outgoing_pattern", 2)
output_incoming <- get_output_path(snakemake@output, "incoming_pattern", 3)
output_signaling_heatmap <- get_output_path(snakemake@output, "signaling_role_heatmap", 4)
output_centrality <- get_output_path(snakemake@output, "centrality_scores", 5)

n_patterns_outgoing <- snakemake@params[["n_patterns_outgoing"]]
n_patterns_incoming <- snakemake@params[["n_patterns_incoming"]]
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

analyze_patterns_single <- function(cellchat, group_name, output_dir,
                                     n_patterns_outgoing, n_patterns_incoming) {
  log_msg(paste("=== Analysis group:", group_name, "==="))
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  pathways <- cellchat@netP$pathways
  n_pathways <- length(pathways)
  n_celltypes <- length(levels(cellchat@idents))
  
  log_msg(paste("  Detected", n_pathways, "signaling pathways"))
  log_msg(paste("  Detected", n_celltypes, "cell types"))
  
  result <- list(
    group = group_name,
    n_pathways = n_pathways,
    n_celltypes = n_celltypes,
    centrality_scores = data.frame()
  )
  
  if (n_pathways < 3) {
    log_msg(paste("  Warning: Insufficient number of pathways (<3), skipping pattern analysis"), level = "WARN")
    return(result)
  }
  
  log_msg("  Computing network centrality scores...")
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  centrality_scores <- data.frame()
  
  if (!is.null(cellchat@netP$centr)) {
    centr_list <- cellchat@netP$centr
    
    for (pathway in names(centr_list)) {
      tryCatch({
        cent <- centr_list[[pathway]]
        if (!is.null(cent) && is.list(cent)) {
          cell_types <- NULL
          for (metric in c("outdeg", "indeg", "flowbet", "info")) {
            if (!is.null(cent[[metric]]) && !is.null(names(cent[[metric]]))) {
              cell_types <- names(cent[[metric]])
              break
            }
          }
          
          if (!is.null(cell_types)) {
            df <- data.frame(
              group = group_name,
              pathway = pathway,
              celltype = cell_types,
              stringsAsFactors = FALSE
            )
            if (!is.null(cent$outdeg)) df$outdeg <- round(as.numeric(cent$outdeg[cell_types]), 6)
            if (!is.null(cent$indeg)) df$indeg <- round(as.numeric(cent$indeg[cell_types]), 6)
            if (!is.null(cent$flowbet)) df$flowbet <- round(as.numeric(cent$flowbet[cell_types]), 6)
            if (!is.null(cent$info)) df$info <- round(as.numeric(cent$info[cell_types]), 6)
            centrality_scores <- rbind(centrality_scores, df)
          }
        }
      }, error = function(e) {})
    }
  }
  
  result$centrality_scores <- centrality_scores
  log_msg(paste("  Extracted", nrow(centrality_scores), "centrality records"))
  
  pdf(file.path(output_dir, paste0("signaling_role_heatmap_", group_name, ".pdf")), 
      width = 14, height = 8)
  
  tryCatch({
    ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",
                                             font.size = 8,
                                             title = paste("Outgoing -", group_name))
    ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",
                                             font.size = 8,
                                             title = paste("Incoming -", group_name))
    print(ht1 + ht2)
    log_msg("  Signaling role heatmap generated successfully")
  }, error = function(e) {
    log_msg(paste("  Failed to generate signaling role heatmap:", e$message), level = "WARN")
    plot.new()
    text(0.5, 0.5, "Failed to generate signaling role heatmap", cex = 1.5)
  })
  
  dev.off()
  
  max_patterns_out <- min(n_patterns_outgoing, n_celltypes - 1, floor(n_pathways / 2))
  max_patterns_out <- max(2, max_patterns_out)
  log_msg(paste("  Outgoing patterns:", max_patterns_out))
  
  pdf(file.path(output_dir, paste0("outgoing_pattern_", group_name, ".pdf")), 
      width = 14, height = 10)
  
  tryCatch({
    cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = max_patterns_out)
    log_msg("  Outgoing communication patterns identified successfully")
    
    tryCatch({
      p_river <- netAnalysis_river(cellchat, pattern = "outgoing")
      print(p_river + ggtitle(paste("Outgoing Patterns -", group_name)))
    }, error = function(e) {})
    
    tryCatch({
      p_dot <- netAnalysis_dot(cellchat, pattern = "outgoing")
      print(p_dot + ggtitle(paste("Outgoing Patterns Dot -", group_name)))
    }, error = function(e) {})
    
  }, error = function(e) {
    log_msg(paste("  Failed to analyze outgoing patterns:", e$message), level = "WARN")
    plot.new()
    text(0.5, 0.5, paste("Failed to analyze outgoing patterns"), cex = 1)
  })
  
  dev.off()
  
  max_patterns_in <- min(n_patterns_incoming, n_celltypes - 1, floor(n_pathways / 2))
  max_patterns_in <- max(2, max_patterns_in)
  log_msg(paste("  Incoming patterns:", max_patterns_in))
  
  pdf(file.path(output_dir, paste0("incoming_pattern_", group_name, ".pdf")), 
      width = 14, height = 10)
  
  tryCatch({
    cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = max_patterns_in)
    log_msg("  Incoming communication patterns identified successfully")
    
    tryCatch({
      p_river <- netAnalysis_river(cellchat, pattern = "incoming")
      print(p_river + ggtitle(paste("Incoming Patterns -", group_name)))
    }, error = function(e) {})
    
    tryCatch({
      p_dot <- netAnalysis_dot(cellchat, pattern = "incoming")
      print(p_dot + ggtitle(paste("Incoming Patterns Dot -", group_name)))
    }, error = function(e) {})
    
  }, error = function(e) {
    log_msg(paste("  Failed to analyze incoming patterns:", e$message), level = "WARN")
    plot.new()
    text(0.5, 0.5, paste("Failed to analyze incoming patterns"), cex = 1)
  })
  
  dev.off()
  
  result$cellchat <- cellchat
  result$max_patterns_out <- max_patterns_out
  result$max_patterns_in <- max_patterns_in
  
  return(result)
}

tryCatch({
  
  log_msg("========== CellChat Communication mode analysis begins ==========")
  
  validate_path <- function(path, name) {
    if (is.null(path) || length(path) == 0 || path == "") {
      stop(paste("Failed to get output path:", name))
    }
  }
  
  validate_path(output_cellchat_final, "cellchat_final")
  validate_path(output_outgoing, "outgoing_pattern")
  validate_path(output_incoming, "incoming_pattern")
  validate_path(output_signaling_heatmap, "signaling_role_heatmap")
  validate_path(output_centrality, "centrality_scores")
  
  for (output_file in c(output_cellchat_final, output_outgoing, output_centrality)) {
    dir_path <- dirname(output_file)
    if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
  }
  
  cellchat_check <- readRDS(input_cellchat)
  if (is.list(cellchat_check) && isTRUE(cellchat_check$skipped)) {
    msg <- ifelse(is.null(cellchat_check$reason), "CellChat步骤已跳过", cellchat_check$reason)
    saveRDS(cellchat_check, output_cellchat_final)
    write_placeholder_pdf(output_outgoing, msg)
    write_placeholder_pdf(output_incoming, msg)
    write_placeholder_pdf(output_signaling_heatmap, msg)
    write.csv(data.frame(message = msg), output_centrality, row.names = FALSE)
    log_msg(paste("Detected placeholder CellChat object, placeholder mode analysis results have been output:", msg), level = "WARN")
    quit(save = "no", status = 0)
  }
  
  cellchat_list_file <- file.path(dirname(input_cellchat), "cellchat_list.rds")
  has_groups <- file.exists(cellchat_list_file)
  
  all_centrality_scores <- data.frame()
  results_list <- list()
  
  if (has_groups) {
    log_msg("Detected grouped CellChat list, performing group-wise analysis...")
    
    cellchat_list <- readRDS(cellchat_list_file)
    group_names <- names(cellchat_list)
    log_msg(paste("Groups:", paste(group_names, collapse = ", ")))
    
    output_base_dir <- dirname(output_outgoing)
    
    for (grp in group_names) {
      group_output_dir <- file.path(output_base_dir, grp)
      
      result <- analyze_patterns_single(
        cellchat = cellchat_list[[grp]],
        group_name = grp,
        output_dir = group_output_dir,
        n_patterns_outgoing = n_patterns_outgoing,
        n_patterns_incoming = n_patterns_incoming
      )
      
      results_list[[grp]] <- result
      
      if (nrow(result$centrality_scores) > 0) {
        all_centrality_scores <- rbind(all_centrality_scores, result$centrality_scores)
      }
      
      if (!is.null(result$cellchat)) {
        cellchat_list[[grp]] <- result$cellchat
      }
    }
    
    if (length(group_names) >= 2) {
      log_msg("\n conduct grouped comparative analysis...")
      
      comparison_dir <- file.path(output_base_dir, "comparison")
      if (!dir.exists(comparison_dir)) dir.create(comparison_dir, recursive = TRUE)
      
      cellchat_merged <- readRDS(input_cellchat)
      
      pdf(file.path(comparison_dir, "signaling_role_comparison.pdf"), 
          width = 16, height = 10)
      
      for (grp in group_names) {
        tryCatch({
          if (!is.null(results_list[[grp]]$cellchat)) {
            cc <- results_list[[grp]]$cellchat
            
            ht1 <- netAnalysis_signalingRole_heatmap(cc, pattern = "outgoing",
                                                     title = paste("Outgoing -", grp))
            print(ht1)
            
            ht2 <- netAnalysis_signalingRole_heatmap(cc, pattern = "incoming",
                                                     title = paste("Incoming -", grp))
            print(ht2)
          }
        }, error = function(e) {
          log_msg(paste("Group", grp, "signaling role heatmap failed:", e$message), level = "WARN")
        })
      }
      
      dev.off()
      
      if (nrow(all_centrality_scores) > 0) {
        pdf(file.path(comparison_dir, "centrality_comparison.pdf"), 
            width = 14, height = 10)
        
        centrality_summary <- all_centrality_scores %>%
          group_by(group, celltype) %>%
          summarise(
            mean_outdeg = mean(outdeg, na.rm = TRUE),
            mean_indeg = mean(indeg, na.rm = TRUE),
            .groups = "drop"
          )
        
        p1 <- ggplot(centrality_summary, aes(x = celltype, y = mean_outdeg, fill = group)) +
          geom_bar(stat = "identity", position = "dodge") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(title = "Mean Outgoing Degree by Cell Type",
               x = "Cell Type", y = "Mean Outgoing Degree")
        print(p1)
        
        p2 <- ggplot(centrality_summary, aes(x = celltype, y = mean_indeg, fill = group)) +
          geom_bar(stat = "identity", position = "dodge") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(title = "Mean Incoming Degree by Cell Type",
               x = "Cell Type", y = "Mean Incoming Degree")
        print(p2)
        
        dev.off()
      }
      
      saveRDS(cellchat_merged, output_cellchat_final)
      log_msg("Merged CellChat object saved")
      
      saveRDS(cellchat_list, file.path(dirname(input_cellchat), "cellchat_list_final.rds"))
      
    } else {
      saveRDS(results_list[[1]]$cellchat, output_cellchat_final)
    }
    
    pdf(output_outgoing, width = 14, height = 10)
    for (grp in group_names) {
      if (!is.null(results_list[[grp]]$cellchat)) {
        tryCatch({
          p <- netAnalysis_river(results_list[[grp]]$cellchat, pattern = "outgoing")
          print(p + ggtitle(paste("Outgoing Patterns -", grp)))
        }, error = function(e) {})
      }
    }
    dev.off()
    
    pdf(output_incoming, width = 14, height = 10)
    for (grp in group_names) {
      if (!is.null(results_list[[grp]]$cellchat)) {
        tryCatch({
          p <- netAnalysis_river(results_list[[grp]]$cellchat, pattern = "incoming")
          print(p + ggtitle(paste("Incoming Patterns -", grp)))
        }, error = function(e) {})
      }
    }
    dev.off()
    
    pdf(output_signaling_heatmap, width = 16, height = 10)
    for (grp in group_names) {
      if (!is.null(results_list[[grp]]$cellchat)) {
        tryCatch({
          ht1 <- netAnalysis_signalingRole_heatmap(results_list[[grp]]$cellchat, 
                                                   pattern = "outgoing",
                                                   title = paste("Outgoing -", grp))
          ht2 <- netAnalysis_signalingRole_heatmap(results_list[[grp]]$cellchat, 
                                                   pattern = "incoming",
                                                   title = paste("Incoming -", grp))
          print(ht1 + ht2)
        }, error = function(e) {})
      }
    }
    dev.off()
    
  } else {
    log_msg("Performing single-group pattern analysis...")
    
    cellchat <- readRDS(input_cellchat)
    
    result <- analyze_patterns_single(
      cellchat = cellchat,
      group_name = "all",
      output_dir = dirname(output_outgoing),
      n_patterns_outgoing = n_patterns_outgoing,
      n_patterns_incoming = n_patterns_incoming
    )
    
    all_centrality_scores <- result$centrality_scores
    
    if (!is.null(result$cellchat)) {
      saveRDS(result$cellchat, output_cellchat_final)
    } else {
      saveRDS(cellchat, output_cellchat_final)
    }
    
    file.copy(file.path(dirname(output_outgoing), "outgoing_pattern_all.pdf"),
              output_outgoing, overwrite = TRUE)
    file.copy(file.path(dirname(output_outgoing), "incoming_pattern_all.pdf"),
              output_incoming, overwrite = TRUE)
    file.copy(file.path(dirname(output_outgoing), "signaling_role_heatmap_all.pdf"),
              output_signaling_heatmap, overwrite = TRUE)
  }
  
  if (nrow(all_centrality_scores) > 0) {
    write.csv(all_centrality_scores, output_centrality, row.names = FALSE, quote = FALSE)
    log_msg(paste("Centrality scores saved (total:", nrow(all_centrality_scores), "records)"))
  } else {
    write.csv(data.frame(message = "Failed to calculate centrality scores"), 
              output_centrality, row.names = FALSE)
    log_msg("Warning: Failed to extract centrality scores", level = "WARN")
  }
  
  cat("\n==========================================\n")
  cat("       CellChat Pattern Analysis Summary\n")
  cat("==========================================\n")
  if (has_groups) {
    cat(paste("Analysis Mode: Group-wise Analysis\n"))
    cat(paste("Number of Groups:", length(group_names), "\n"))
    cat(paste("Groups:", paste(group_names, collapse = ", "), "\n"))
  } else {
    cat("Analysis Mode: Overall Analysis\n")
  }
  cat("==========================================\n\n")
  
  log_msg("========== CellChat Pattern Analysis Complete ==========")
  
}, error = function(e) {
  log_msg(paste("Error:", e$message), level = "ERROR")
  log_msg(paste("Error Traceback:", paste(capture.output(traceback()), collapse = "\n")), level = "ERROR")
  stop(e)
})
