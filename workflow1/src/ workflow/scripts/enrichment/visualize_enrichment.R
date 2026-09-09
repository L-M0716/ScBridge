suppressPackageStartupMessages({
  library(ggplot2)
  library(clusterProfiler)
  library(dplyr)
  library(patchwork)
  library(enrichplot)
  library(tidyr)
})


rds_file <- snakemake@input[["enrichment_rds"]]

top_n <- snakemake@params[["top_n_terms"]]
plot_width <- snakemake@params[["plot_width"]]
plot_height <- snakemake@params[["plot_height"]]

go_viz_dir <- snakemake@output[["go_viz_dir"]]
go_done <- snakemake@output[["go_done"]]
kegg_viz_dir <- snakemake@output[["kegg_viz_dir"]]
kegg_done <- snakemake@output[["kegg_done"]]

for (dir in c(go_viz_dir, kegg_viz_dir)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}


create_enhanced_dotplot <- function(enrich_result, title, n_terms = top_n) {
  if (is.null(enrich_result) || nrow(enrich_result) == 0) {
    return(NULL)
  }
  
  n_terms <- min(n_terms, nrow(enrich_result))
  
  p <- dotplot(enrich_result, showCategory = n_terms) +
    ggtitle(title) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9),
      axis.title = element_text(size = 10),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    ) +
    scale_color_gradient(low = "red", high = "blue", name = "Adj.P")
  
  return(p)
}

create_enhanced_barplot <- function(enrich_result, title, n_terms = top_n) {
  if (is.null(enrich_result) || nrow(enrich_result) == 0) {
    return(NULL)
  }
  
  n_terms <- min(n_terms, nrow(enrich_result))
  plot_data <- as.data.frame(enrich_result)[1:n_terms, ]
  
  plot_data$significance <- cut(plot_data$p.adjust, 
                                breaks = c(0, 0.001, 0.01, 0.05, 1),
                                labels = c("***", "**", "*", "ns"))
  
  plot_data$neg_log_fdr <- -log10(plot_data$p.adjust)
  
  plot_data$Description_short <- ifelse(
    nchar(plot_data$Description) > 50,
    paste0(substr(plot_data$Description, 1, 47), "..."),
    plot_data$Description
  )
  
  p <- ggplot(plot_data, aes(x = Count, y = reorder(Description_short, Count))) +
    geom_bar(stat = "identity", aes(fill = neg_log_fdr)) +
    geom_text(aes(label = significance), hjust = -0.2, size = 3) +
    scale_fill_gradient(low = "lightblue", high = "darkred", 
                        name = "-log10(FDR)") +
    labs(title = title, x = "Gene Count", y = "") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 9),
      axis.title = element_text(size = 10),
      legend.position = "right",
      panel.grid.major.y = element_blank()
    ) +
    xlim(0, max(plot_data$Count) * 1.15)
  
  return(p)
}


create_celltype_comparison_heatmap <- function(all_results, comparison_filter, 
                                                ontology = "BP", top_n = 10) {
  combined_data <- data.frame()
  
  for (result_key in names(all_results)) {
    result_info <- all_results[[result_key]]
    comparison <- result_info$comparison
    cell_type <- result_info$cell_type
    
    if (comparison != comparison_filter) next
    
    if (!"all" %in% names(result_info$results)) next
    
    result <- result_info$results[["all"]]
    
    if (ontology == "KEGG") {
      enrich_res <- result$kegg
    } else {
      enrich_res <- result$go[[ontology]]
    }
    
    if (!is.null(enrich_res) && nrow(enrich_res) > 0) {
      df <- as.data.frame(enrich_res)[1:min(top_n, nrow(enrich_res)), ]
      df$cell_type <- cell_type
      df$neg_log_p <- -log10(df$p.adjust)
      combined_data <- rbind(combined_data, df[, c("Description", "cell_type", "neg_log_p", "Count")])
    }
  }
  
  if (nrow(combined_data) == 0) return(NULL)
  
  top_terms <- combined_data %>%
    group_by(Description) %>%
    summarise(n = n(), mean_sig = mean(neg_log_p)) %>%
    arrange(desc(n), desc(mean_sig)) %>%
    head(20) %>%
    pull(Description)
  
  plot_data <- combined_data %>%
    filter(Description %in% top_terms) %>%
    complete(Description, cell_type, fill = list(neg_log_p = 0, Count = 0))
  
  plot_data$Description_short <- ifelse(
    nchar(plot_data$Description) > 40,
    paste0(substr(plot_data$Description, 1, 37), "..."),
    plot_data$Description
  )
  
  p <- ggplot(plot_data, aes(x = cell_type, y = Description_short, fill = neg_log_p)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(Count > 0, Count, "")), size = 3) +
    scale_fill_gradient(low = "white", high = "red", name = "-log10(FDR)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 9),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(title = paste(ontology, "Enrichment Across Cell Types\n(", comparison_filter, ")"))
  
  return(p)
}
create_updown_comparison_plot <- function(all_results, comparison_filter, 
                                          cell_type_filter, ontology = "BP", top_n = 10) {
  combined_data <- data.frame()
  
  for (result_key in names(all_results)) {
    result_info <- all_results[[result_key]]
    comparison <- result_info$comparison
    cell_type <- result_info$cell_type
    
    if (comparison != comparison_filter || cell_type != cell_type_filter) next
    
    for (direction in c("up", "down")) {
      if (!direction %in% names(result_info$results)) next
      
      result <- result_info$results[[direction]]
      
      if (ontology == "KEGG") {
        enrich_res <- result$kegg
      } else {
        enrich_res <- result$go[[ontology]]
      }
      
      if (!is.null(enrich_res) && nrow(enrich_res) > 0) {
        df <- as.data.frame(enrich_res)[1:min(top_n, nrow(enrich_res)), ]
        df$direction <- direction
        df$neg_log_p <- -log10(df$p.adjust)
        if (direction == "down") {
          df$neg_log_p <- -df$neg_log_p
        }
        combined_data <- rbind(combined_data, df[, c("Description", "direction", "neg_log_p", "Count")])
      }
    }
  }
  
  if (nrow(combined_data) == 0) return(NULL)
  
  combined_data$Description_short <- ifelse(
    nchar(combined_data$Description) > 45,
    paste0(substr(combined_data$Description, 1, 42), "..."),
    combined_data$Description
  )
  
  p <- ggplot(combined_data, aes(x = neg_log_p, y = reorder(Description_short, neg_log_p), 
                                  fill = direction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"),
                      labels = c("up" = "Up-regulated", "down" = "Down-regulated")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "top"
    ) +
    labs(title = paste(ontology, "-", cell_type_filter, "\n(", comparison_filter, ")"),
         x = "-log10(FDR) * direction",
         y = "",
         fill = "Regulation")
  
  return(p)
}
tryCatch({
  message("==========================================================")
  message("Start pathway enrichment result visualization")
  message("==========================================================")
  
  if (!file.exists(rds_file)) {
    stop(paste("Enrichment analysis result file not found:", rds_file))
  }
  
  all_results <- readRDS(rds_file)
  
  all_comparisons <- unique(sapply(all_results, function(x) x$comparison))
  all_cell_types <- unique(sapply(all_results, function(x) x$cell_type))
  
  message(paste("Found", length(all_comparisons), "comparison groups:", 
                paste(all_comparisons, collapse = ", ")))
  message(paste("Found", length(all_cell_types), "cell types:", 
                paste(all_cell_types, collapse = ", ")))
  
  plot_count <- 0
  for (comparison in all_comparisons) {
    comparison_safe <- gsub("[^A-Za-z0-9_-]", "_", comparison)
    
    comp_go_viz_dir <- file.path(go_viz_dir, comparison_safe)
    comp_kegg_viz_dir <- file.path(kegg_viz_dir, comparison_safe)
    dir.create(comp_go_viz_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(comp_kegg_viz_dir, recursive = TRUE, showWarnings = FALSE)
    
    message(paste("\nProcessing comparison group:", comparison))
    

    for (cell_type in all_cell_types) {
      cell_type_safe <- gsub("[^A-Za-z0-9_-]", "_", cell_type)
      result_key <- paste(comparison_safe, cell_type_safe, sep = "__")
      
      if (!result_key %in% names(all_results)) next
      
      result_info <- all_results[[result_key]]
      cell_results <- result_info$results
      
      for (gene_list in names(cell_results)) {
        result <- cell_results[[gene_list]]
        
        for (ontology in c("BP", "MF", "CC")) {
          if (!is.null(result$go[[ontology]]) && nrow(result$go[[ontology]]) > 0) {
            dot_file <- file.path(comp_go_viz_dir, 
                                  paste0(cell_type_safe, "_", gene_list, "_GO_", ontology, "_dot.pdf"))
            bar_file <- file.path(comp_go_viz_dir, 
                                  paste0(cell_type_safe, "_", gene_list, "_GO_", ontology, "_bar.pdf"))
            
            title <- paste(cell_type, "-", gene_list, "\nGO", ontology)
            
            p_dot <- create_enhanced_dotplot(result$go[[ontology]], title)
            p_bar <- create_enhanced_barplot(result$go[[ontology]], title)
            
            if (!is.null(p_dot)) {
              ggsave(dot_file, p_dot, width = plot_width, height = plot_height, 
                     device = "pdf", useDingbats = FALSE)
              plot_count <- plot_count + 1
            }
            if (!is.null(p_bar)) {
              ggsave(bar_file, p_bar, width = plot_width, height = plot_height, 
                     device = "pdf", useDingbats = FALSE)
              plot_count <- plot_count + 1
            }
          }
        }
        
        if (!is.null(result$kegg) && nrow(result$kegg) > 0) {
          kegg_dot_file <- file.path(comp_kegg_viz_dir, 
                                     paste0(cell_type_safe, "_", gene_list, "_KEGG_dot.pdf"))
          kegg_bar_file <- file.path(comp_kegg_viz_dir, 
                                     paste0(cell_type_safe, "_", gene_list, "_KEGG_bar.pdf"))
          
          title <- paste(cell_type, "-", gene_list, "\nKEGG")
          
          p_kegg_dot <- create_enhanced_dotplot(result$kegg, title)
          p_kegg_bar <- create_enhanced_barplot(result$kegg, title)
          
          if (!is.null(p_kegg_dot)) {
            ggsave(kegg_dot_file, p_kegg_dot, width = plot_width, height = plot_height, 
                   device = "pdf", useDingbats = FALSE)
            plot_count <- plot_count + 1
          }
          if (!is.null(p_kegg_bar)) {
            ggsave(kegg_bar_file, p_kegg_bar, width = plot_width, height = plot_height, 
                   device = "pdf", useDingbats = FALSE)
            plot_count <- plot_count + 1
          }
        }
      }
      

      for (ontology in c("BP", "KEGG")) {
        p_updown <- create_updown_comparison_plot(all_results, comparison, cell_type, ontology)
        
        if (!is.null(p_updown)) {
          if (ontology == "BP") {
            updown_file <- file.path(comp_go_viz_dir, 
                                     paste0(cell_type_safe, "_updown_comparison_GO_BP.pdf"))
          } else {
            updown_file <- file.path(comp_kegg_viz_dir, 
                                     paste0(cell_type_safe, "_updown_comparison_KEGG.pdf"))
          }
          
          ggsave(updown_file, p_updown, width = 10, height = 8, 
                 device = "pdf", useDingbats = FALSE)
          plot_count <- plot_count + 1
          message(paste("  Saving up/down comparison plot:", basename(updown_file)))
        }
      }
    }
    

    message(paste("  Creating cross-cell type comparison heatmap..."))
    
    for (ontology in c("BP", "KEGG")) {
      p_heatmap <- create_celltype_comparison_heatmap(all_results, comparison, ontology)
      
      if (!is.null(p_heatmap)) {
        if (ontology == "BP") {
          heatmap_file <- file.path(comp_go_viz_dir, "celltype_comparison_heatmap_GO_BP.pdf")
        } else {
          heatmap_file <- file.path(comp_kegg_viz_dir, "celltype_comparison_heatmap_KEGG.pdf")
        }
        
        ggsave(heatmap_file, p_heatmap, width = 12, height = 10, 
               device = "pdf", useDingbats = FALSE)
        plot_count <- plot_count + 1
        message(paste("  Saving cross-cell type heatmap:", basename(heatmap_file)))
      }
    }
  }
  

  if (length(all_comparisons) > 1) {
    message("\nCreating global comprehensive comparison plot...")
    

    global_data <- data.frame()
    
    for (result_key in names(all_results)) {
      result_info <- all_results[[result_key]]
      if (!"all" %in% names(result_info$results)) next
      
      result <- result_info$results[["all"]]
      if (!is.null(result$go$BP) && nrow(result$go$BP) > 0) {
        global_data <- rbind(global_data, data.frame(
          comparison = result_info$comparison,
          cell_type = result_info$cell_type,
          n_terms = nrow(result$go$BP),
          top_term = as.data.frame(result$go$BP)$Description[1]
        ))
      }
    }
    
    if (nrow(global_data) > 0) {
      p_global <- ggplot(global_data, aes(x = comparison, y = cell_type, fill = n_terms)) +
        geom_tile(color = "white") +
        geom_text(aes(label = n_terms), size = 4) +
        scale_fill_gradient(low = "lightyellow", high = "darkred", name = "# GO BP Terms") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_blank()
        ) +
        labs(title = "Number of Significant GO BP Terms\nAcross Comparisons and Cell Types")
      
      global_file <- file.path(go_viz_dir, "global_summary_heatmap.pdf")
      ggsave(global_file, p_global, width = 10, height = 8, 
             device = "pdf", useDingbats = FALSE)
      plot_count <- plot_count + 1
      message(paste("Saving global comprehensive heatmap:", basename(global_file)))
    }
  }
  

  message("\n==========================================================")
  message(paste("Pathway enrichment visualization completed!"))
  message(paste("Total PDF plots generated:", plot_count))
  message(paste("GO visualization results located at:", go_viz_dir))
  message(paste("KEGG visualization results located at:", kegg_viz_dir))
  message("==========================================================")
  
  writeLines(paste("Visualization completed at", Sys.time()), go_done)
  writeLines(paste("Visualization completed at", Sys.time()), kegg_done)
  
}, error = function(e) {
  message(paste("Error:", e$message))
  print(traceback())
  stop(e)
})