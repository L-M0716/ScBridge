suppressPackageStartupMessages({
  library(clusterProfiler)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(DOSE)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
})

de_results_file <- snakemake@input[["de_results"]]
seurat_file <- snakemake@input[["final_annotated_seurat"]]

go_dir   <- snakemake@output[["go_dir"]]
kegg_dir <- snakemake@output[["kegg_dir"]]
base_output_dir <- dirname(snakemake@output[["enrichment_summary"]])
for (dir in c(go_dir, kegg_dir)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}
output_files <- list(
  rds_all = snakemake@output[["enrichment_rds"]],
  summary = snakemake@output[["enrichment_summary"]],
  summary_text = snakemake@output[["enrichment_summary_txt"]],
  unmapped_genes = snakemake@output[["unmapped_genes"]],
  log = snakemake@output[["enrichment_log"]]
)

pval_cutoff <- as.numeric(snakemake@params[["pval_cutoff"]])
qval_cutoff <- as.numeric(snakemake@params[["qval_cutoff"]])
logfc_threshold <- as.numeric(snakemake@params[["logfc_threshold"]])
min_gene_size <- as.integer(snakemake@params[["min_gene_size"]])
max_gene_size <- as.integer(snakemake@params[["max_gene_size"]])
organism <- snakemake@params[["organism"]]
human_kegg_db <- snakemake@params[["human_db_path"]]
mouse_kegg_db <- snakemake@params[["mouse_db_path"]]


log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] %s", timestamp, level, msg)
  cat(line, "\n", file = output_files$log, append = TRUE)
  message(line)
}
setup_organism_databases <- function(organism) {
  organism_lower <- tolower(organism)
  
  if (organism_lower %in% c("mouse", "mm", "mus musculus", "mus_musculus", "mmu")) {
    log_msg("Species detected: Mus musculus (Mouse)")
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    return(list(
      org_db = org.Mm.eg.db,
      kegg_organism = "mmu",
      kegg_db_path = mouse_kegg_db,
      species_name = "mouse"
    ))
  } else if (organism_lower %in% c("human", "hs", "homo sapiens", "homo_sapiens", "hsa")) {
    log_msg("Species detected: Homo sapiens (Human)")
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    return(list(
      org_db = org.Hs.eg.db,
      kegg_organism = "hsa",
      kegg_db_path = human_kegg_db,
      species_name = "human"
    ))
  } else {
    stop(paste("Unsupported species:", organism, " please use 'mouse' or 'human'"))
  }
}
organism_info <- setup_organism_databases(organism)
org_db <- organism_info$org_db
kegg_organism <- organism_info$kegg_organism
kegg_offline_data_file <- organism_info$kegg_db_path

log_msg(paste("Using database:", organism_info$species_name))
log_msg(paste("KEGG database path:", kegg_offline_data_file))


load_offline_kegg_data <- function() {
  if (is.null(kegg_offline_data_file) || !file.exists(kegg_offline_data_file)) {
    log_msg(paste("Warning: Offline KEGG data file does not exist:", kegg_offline_data_file), level = "WARN")
    return(NULL)
  }
  
  log_msg(paste("Loading offline KEGG data:", kegg_offline_data_file))
  kegg_data <- readRDS(kegg_offline_data_file)
  
  log_msg(paste("KEGG data loaded successfully - Pathways:", length(unique(kegg_data$term2gene$term)),
                "Genes:", length(unique(kegg_data$term2gene$gene))))
  
  return(kegg_data)
}


convert_gene_symbols <- function(gene_symbols, context = "") {
  gene_symbols <- unique(na.omit(gene_symbols))
  
  if (length(gene_symbols) == 0) {
    log_msg(paste("Warning: No valid gene symbols found", context), level = "WARN")
    return(character(0))
  }
  
  converted <- tryCatch({
    bitr(gene_symbols, 
         fromType = "SYMBOL", 
         toType = "ENTREZID", 
         OrgDb = org_db)
  }, error = function(e) {
    log_msg(paste("Warning: Gene ID conversion failed:", e$message), level = "ERROR")
    return(NULL)
  })
  
  if (is.null(converted) || nrow(converted) == 0) return(character(0))
  
  unmapped <- setdiff(gene_symbols, converted$SYMBOL)
  if (length(unmapped) > 0) {
    log_msg(paste("Warning:", context, "中", length(unmapped), 
                  "genes cannot be mapped to an Entrez ID"), level = "WARN")
    
    unmapped_df <- data.frame(
      gene_symbol = unmapped,
      context = context,
      organism = organism_info$species_name,
      timestamp = Sys.time()
    )
    
    if (file.exists(output_files$unmapped_genes)) {
      write.table(unmapped_df, file = output_files$unmapped_genes, 
                  append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
    } else {
      write.csv(unmapped_df, file = output_files$unmapped_genes, row.names = FALSE)
    }
  }
  
  entrez_ids <- unique(converted$ENTREZID)
  log_msg(paste("Successfully converted", length(entrez_ids), "unique Entrez IDs", context))
  
  return(entrez_ids)
}


run_go_enrichment <- function(entrez_ids, ont = "BP", context = "") {
  if (length(entrez_ids) < 5) {
    log_msg(paste("Not enough genes (<5), skipping GO", ont, "analysis", context), level = "WARN")
    return(NULL)
  }
  
  log_msg(paste("Executing GO", ont, "enrichment analysis", context, "- Number of genes:", length(entrez_ids)))
  
  go_result <- tryCatch({
    enrichGO(
      gene = entrez_ids,
      OrgDb = org_db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = pval_cutoff,
      qvalueCutoff = qval_cutoff,
      minGSSize = min_gene_size,
      maxGSSize = max_gene_size,
      readable = TRUE
    )
  }, error = function(e) {
    log_msg(paste("GO enrichment analysis failed:", e$message), level = "ERROR")
    return(NULL)
  })
  
  if (!is.null(go_result) && nrow(go_result) > 0) {
    log_msg(paste("GO", ont, "found", nrow(go_result), "significant pathways"))
  }
  
  return(go_result)
}


run_kegg_enrichment_offline <- function(entrez_ids, kegg_data, context = "") {
  if (length(entrez_ids) < 5) {
    log_msg(paste("Not enough genes (<5), skipping KEGG analysis", context), level = "WARN")
    return(NULL)
  }
  
  if (is.null(kegg_data)) {
    log_msg("Offline KEGG data is not available, skipping KEGG analysis", level = "WARN")
    return(NULL)
  }
  
  log_msg(paste("Executing offline KEGG pathway enrichment analysis", context, "- Number of genes:", length(entrez_ids)))
  
  kegg_result <- tryCatch({
    enricher(
      gene = entrez_ids,
      TERM2GENE = kegg_data$term2gene,
      TERM2NAME = kegg_data$term2name,
      pAdjustMethod = "BH",
      pvalueCutoff = pval_cutoff,
      qvalueCutoff = qval_cutoff,
      minGSSize = min_gene_size,
      maxGSSize = max_gene_size
    )
  }, error = function(e) {
    log_msg(paste("KEGG enrichment analysis failed:", e$message), level = "ERROR")
    return(NULL)
  })
  
  if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
    log_msg(paste("KEGG found", nrow(kegg_result), "significant pathways"))
    
    kegg_result <- tryCatch({
      all_genes <- unique(unlist(strsplit(kegg_result@result$geneID, "/")))
      
      gene_mapping <- bitr(all_genes,
                           fromType = "ENTREZID",
                           toType = "SYMBOL",
                           OrgDb = org_db)
      
      id_to_symbol <- setNames(gene_mapping$SYMBOL, gene_mapping$ENTREZID)
      
      kegg_result@result$geneID <- sapply(kegg_result@result$geneID, function(x) {
        ids <- unlist(strsplit(x, "/"))
        symbols <- id_to_symbol[ids]
        paste(na.omit(symbols), collapse = "/")
      })
      
      return(kegg_result)
    }, error = function(e) {
      log_msg("KEGG result conversion to readable format failed, keeping Entrez IDs", level = "WARN")
      return(kegg_result)
    })
  }
  
  return(kegg_result)
}


save_enrichment_results <- function(result, file_path, result_type) {
  if (!is.null(result) && nrow(result) > 0) {
    tryCatch({
      dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
      write.csv(as.data.frame(result), file = file_path, row.names = FALSE)
      log_msg(paste("save", result_type, "results:", nrow(result), "entries ->", basename(file_path)))
      return(TRUE)
    }, error = function(e) {
      log_msg(paste("Failed to save results:", e$message), level = "ERROR")
      return(FALSE)
    })
  } else {
    log_msg(paste("No significant", result_type, "results found"), level = "WARN")
    return(FALSE)
  }
}


get_group_info_from_seurat <- function(seurat_file) {
  log_msg(paste("Reading Seurat object:", seurat_file))
  
  if (!file.exists(seurat_file)) {
    stop(paste("Seurat file does not exist:", seurat_file))
  }
  
  seurat_obj <- readRDS(seurat_file)
  meta_data <- seurat_obj@meta.data
  
  group_info <- list()
  
  if ("group" %in% colnames(meta_data)) {
    group_info$groups <- unique(meta_data$group)
    log_msg(paste("Found group column, groups:", paste(group_info$groups, collapse = ", ")))
  } else {
    log_msg("Warning: group column not found in Seurat object", level = "WARN")
    group_info$groups <- NULL
  }
  
  cell_type_cols <- c("cell_type", "celltype", "CellType", "cell_type_annotation", 
                      "seurat_clusters", "cluster", "ident")
  
  for (col in cell_type_cols) {
    if (col %in% colnames(meta_data)) {
      group_info$cell_types <- unique(meta_data[[col]])
      group_info$cell_type_col <- col
      log_msg(paste("Found cell type column:", col, "- Number of types:", length(group_info$cell_types)))
      break
    }
  }
  
  if (is.null(group_info$cell_type_col)) {
    log_msg("Warning: Cell type annotation column not found", level = "WARN")
  }
  
  rm(seurat_obj)
  gc()
  
  return(group_info)
}


tryCatch({
  log_msg("==========================================================")
  log_msg("Starting pathway enrichment analysis workflow")
  log_msg("==========================================================")
  log_msg(paste("Species:", organism_info$species_name))
  
  kegg_data <- load_offline_kegg_data()
  
  group_info <- get_group_info_from_seurat(seurat_file)
  
  log_msg(paste("Loading differential expression results:", de_results_file))
  
  if (!file.exists(de_results_file)) {
    stop("Differential expression results file does not exist")
  }
  
  de_results <- read.csv(de_results_file)
  
  if (nrow(de_results) == 0) {
    stop("Differential expression results are empty")
  }
  
  log_msg(paste("Differential expression results: ", nrow(de_results), "rows,", ncol(de_results), "columns"))
  log_msg(paste("Column names:", paste(colnames(de_results), collapse = ", ")))
  
  required_cols <- c("gene", "p_val_adj", "avg_log2FC")
  missing_cols <- setdiff(required_cols, colnames(de_results))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  

  
  comparison_col <- NULL
  if ("comparison" %in% colnames(de_results)) {
    comparison_col <- "comparison"
  } else if ("contrast" %in% colnames(de_results)) {
    comparison_col <- "contrast"
  } else if ("group" %in% colnames(de_results)) {
    comparison_col <- "group"
  }
  
  if (!is.null(comparison_col)) {
    comparisons <- unique(de_results[[comparison_col]])
    de_results$comparison <- de_results[[comparison_col]]
    log_msg(paste("Discovery Comparison Series:", comparison_col, "- Comparisons:", paste(comparisons, collapse = ", ")))
  } else {
    if (!is.null(group_info$groups) && length(group_info$groups) == 2) {
      comparison_name <- paste(group_info$groups[1], "vs", group_info$groups[2], sep = "_")
      de_results$comparison <- comparison_name
      comparisons <- comparison_name
      log_msg(paste("Construct comparison based on Seurat grouping information:", comparison_name))
    } else {
      de_results$comparison <- "DEG_analysis"
      comparisons <- "DEG_analysis"
      log_msg("Warning: Comparison information not found, using default label")
    }
  }

  cell_type_col <- NULL
  for (col in c("cell_type", "celltype", "CellType", "cluster", "seurat_clusters")) {
    if (col %in% colnames(de_results)) {
      cell_type_col <- col
      break
    }
  }
  
  if (!is.null(cell_type_col)) {
    cell_types <- unique(de_results[[cell_type_col]])
    de_results$cell_type <- de_results[[cell_type_col]]
    log_msg(paste("Found cell type column:", cell_type_col, "- Number of types:", length(cell_types)))
  } else {
    de_results$cell_type <- "all_cells"
    cell_types <- "all_cells"
    log_msg("Warning: Cell type annotation column not found, treating as all cells")
  }
  
  log_msg(paste("\nAnalysis strategy: Enrichment analysis by [Comparison x Cell Type]"))
  log_msg(paste("Number of comparisons:", length(comparisons)))
  log_msg(paste("Number of cell types:", length(cell_types)))
  log_msg(paste("Expected analysis groups:", length(comparisons) * length(cell_types) * 3, "(including all/up/down)"))
  

  
  all_enrichment_results <- list()
  analysis_summary <- data.frame()
  
  for (comparison in comparisons) {
    log_msg(paste("\n", paste(rep("=", 60), collapse = "")))
    log_msg(paste("Processing comparison:", comparison))
    log_msg(paste(rep("=", 60), collapse = ""))
    

    comparison_safe <- gsub("[^A-Za-z0-9_-]", "_", comparison)
    comparison_go_dir <- file.path(go_dir, comparison_safe)
    comparison_kegg_dir <- file.path(kegg_dir, comparison_safe)
    dir.create(comparison_go_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(comparison_kegg_dir, recursive = TRUE, showWarnings = FALSE)
    
    for (cell_type in cell_types) {
      log_msg(paste("\n--- Cell Type:", cell_type, "---"))
      
      cell_de <- de_results %>% 
        filter(comparison == !!comparison, cell_type == !!cell_type)
      
      if (nrow(cell_de) == 0) {
        log_msg(paste("Skipping: No differentially expressed gene data"), level = "WARN")
        next
      }
      
      log_msg(paste("This combination has", nrow(cell_de), "gene records"))
      
      all_sig <- cell_de %>% 
        filter(p_val_adj < pval_cutoff, abs(avg_log2FC) > logfc_threshold)
      
      up_genes <- all_sig %>% 
        filter(avg_log2FC > 0) %>% 
        pull(gene)
      
      down_genes <- all_sig %>% 
        filter(avg_log2FC < 0) %>% 
        pull(gene)
      
      all_genes <- all_sig %>% pull(gene)
      
      log_msg(sprintf("Significantly differentially expressed genes: Total=%d, Up=%d, Down=%d",
                      length(all_genes), length(up_genes), length(down_genes)))
      
      gene_lists <- list(
        "all" = all_genes,
        "up" = up_genes,
        "down" = down_genes
      )
      
      cell_type_safe <- gsub("[^A-Za-z0-9_-]", "_", cell_type)
      cell_results <- list()
      
      for (list_name in names(gene_lists)) {
        genes <- gene_lists[[list_name]]
        context <- paste("[", comparison, "|", cell_type, "|", list_name, "]")
        
        if (length(genes) < 5) {
          log_msg(paste("Skipping", list_name, ": Insufficient genes (", length(genes), "<5)"), level = "WARN")
          
          summary_row <- data.frame(
            comparison = comparison,
            cell_type = cell_type,
            gene_list = list_name,
            organism = organism_info$species_name,
            n_genes = length(genes),
            n_entrez = 0,
            go_bp = 0, go_mf = 0, go_cc = 0, kegg = 0,
            status = "skipped_insufficient_genes",
            timestamp = Sys.time()
          )
          analysis_summary <- rbind(analysis_summary, summary_row)
          next
        }
        
        log_msg(paste("Processing", context, ":", length(genes), "genes"))
        
        entrez_ids <- convert_gene_symbols(genes, context)
        
        if (length(entrez_ids) < 5) {
          log_msg(paste("Skipping", context, ": Insufficient valid Entrez IDs"), level = "WARN")
          
          summary_row <- data.frame(
            comparison = comparison,
            cell_type = cell_type,
            gene_list = list_name,
            organism = organism_info$species_name,
            n_genes = length(genes),
            n_entrez = length(entrez_ids),
            go_bp = 0, go_mf = 0, go_cc = 0, kegg = 0,
            status = "skipped_insufficient_entrez",
            timestamp = Sys.time()
          )
          analysis_summary <- rbind(analysis_summary, summary_row)
          next
        }
        
        go_results <- list()
        for (ontology in c("BP", "MF", "CC")) {
          output_file <- file.path(comparison_go_dir, 
                                   paste0(cell_type_safe, "_", list_name, "_GO_", ontology, ".csv"))
          
          go_result <- run_go_enrichment(entrez_ids, ont = ontology, context = context)
          
          if (save_enrichment_results(go_result, output_file, paste("GO", ontology))) {
            go_results[[ontology]] <- go_result
          }
        }
        
        kegg_output <- file.path(comparison_kegg_dir, 
                                 paste0(cell_type_safe, "_", list_name, "_KEGG.csv"))
        kegg_result <- run_kegg_enrichment_offline(entrez_ids, kegg_data, context = context)
        save_enrichment_results(kegg_result, kegg_output, "KEGG")
        
        cell_results[[list_name]] <- list(
          genes = genes,
          entrez_ids = entrez_ids,
          go = go_results,
          kegg = kegg_result
        )
        
        summary_row <- data.frame(
          comparison = comparison,
          cell_type = cell_type,
          gene_list = list_name,
          organism = organism_info$species_name,
          n_genes = length(genes),
          n_entrez = length(entrez_ids),
          go_bp = if(!is.null(go_results$BP)) nrow(go_results$BP) else 0,
          go_mf = if(!is.null(go_results$MF)) nrow(go_results$MF) else 0,
          go_cc = if(!is.null(go_results$CC)) nrow(go_results$CC) else 0,
          kegg = if(!is.null(kegg_result)) nrow(kegg_result) else 0,
          status = "completed",
          timestamp = Sys.time()
        )
        analysis_summary <- rbind(analysis_summary, summary_row)
      }
      
      result_key <- paste(comparison_safe, cell_type_safe, sep = "__")
      all_enrichment_results[[result_key]] <- list(
        comparison = comparison,
        cell_type = cell_type,
        results = cell_results
      )
    }
  }
  

  
  log_msg("\nSave comprehensive results...")
  
  saveRDS(all_enrichment_results, output_files$rds_all)
  log_msg(paste("RDS results saved:", output_files$rds_all))
  
  write.csv(analysis_summary, output_files$summary, row.names = FALSE)
  log_msg(paste("Analysis summary saved:", output_files$summary))
  
  sink(output_files$summary_text)
  cat("============================================================\n")
  cat("Pathway Enrichment Analysis Summary\n")
  cat("Analysis Strategy: By Comparison + Cell Type\n")
  cat("============================================================\n\n")
  
  cat("Analysis Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Input Files:\n")
  cat("  - Differentially Expressed Gene Results:", de_results_file, "\n")
  cat("  - Seurat Object:", seurat_file, "\n\n")
  
  cat("Organism Information:\n")
  cat("  - Species:", organism_info$species_name, "\n")
  cat("  - KEGG Species Code:", kegg_organism, "\n")
  
  if (!is.null(kegg_data)) {
    cat("  - KEGG Data Source: Offline Data\n")
    cat("  - KEGG Data File:", kegg_offline_data_file, "\n")
    if (!is.null(kegg_data$download_date)) {
      cat("  - KEGG Data Date:", as.character(kegg_data$download_date), "\n")
    }
  }
  
  cat("\nParameter Settings:\n")
  cat("  - P-value Threshold:", pval_cutoff, "\n")
  cat("  - Q-value Threshold:", qval_cutoff, "\n")
  cat("  - LogFC Threshold:", logfc_threshold, "\n")
  cat("  - Gene Set Size:", min_gene_size, "-", max_gene_size, "\n")
  
  cat("\nAnalysis Dimensions:\n")
  cat("  - Comparisons:", paste(comparisons, collapse = ", "), "\n")
  cat("  - Cell Types:", paste(cell_types, collapse = ", "), "\n")
  

  
  if (nrow(analysis_summary) > 0) {
    cat("Statistics by comparison group:\n")
    comp_stats <- analysis_summary %>%
      group_by(comparison) %>%
      summarise(
        n_analyses = n(),
        completed = sum(status == "completed"),
        total_go_bp = sum(go_bp),
        total_kegg = sum(kegg)
      )
    print(as.data.frame(comp_stats))
    
    cat("\nStatistics by cell type:\n")
    cell_stats <- analysis_summary %>%
      group_by(cell_type) %>%
      summarise(
        n_analyses = n(),
        completed = sum(status == "completed"),
        total_go_bp = sum(go_bp),
        total_kegg = sum(kegg)
      )
    print(as.data.frame(cell_stats))
    
    cat("\nDetailed Results:\n")
    print(as.data.frame(analysis_summary))
  }
  
  cat("GO Results:", go_dir, "\n")
  for (comp in unique(analysis_summary$comparison)) {
    cat("  └──", gsub("[^A-Za-z0-9_-]", "_", comp), "/\n")
  }
  cat("KEGG Results:", kegg_dir, "\n")
  for (comp in unique(analysis_summary$comparison)) {
    cat("  └──", gsub("[^A-Za-z0-9_-]", "_", comp), "/\n")
  }
  sink()
  

  log_msg(paste("Analysis completed for", length(comparisons), "comparison groups x", 
                length(cell_types), "cell types"))
  log_msg(paste("Completed analyses:", sum(analysis_summary$status == "completed")))
  log_msg(paste("Total significant GO BP pathways:", sum(analysis_summary$go_bp)))
  log_msg(paste("Total significant KEGG pathways:", sum(analysis_summary$kegg)))
  
}, error = function(e) {
  log_msg(paste("Fatal error:", e$message), level = "ERROR")
  log_msg(paste("Error details:", paste(capture.output(traceback()), collapse = "\n")), level = "ERROR")
  stop(e)
})