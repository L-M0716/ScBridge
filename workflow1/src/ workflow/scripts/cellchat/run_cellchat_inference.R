# scripts/cellchat/run_cellchat_inference.R

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(tidyverse)
  library(Matrix)
})

input_file <- snakemake@input[["annotated"]]
output_cellchat <- snakemake@output[["cellchat_obj"]]
output_net_df <- snakemake@output[["net_df"]]
output_lr_pairs <- snakemake@output[["lr_pairs"]]
db_type <- snakemake@params[["db_type"]]
min_cells <- snakemake@params[["min_cells"]]
cell_type_column <- snakemake@params[["cell_type_column"]]
organism <- snakemake@params[["organism"]]

log_file <- snakemake@log[[1]]

if (!dir.exists(dirname(log_file))) dir.create(dirname(log_file), recursive = TRUE)
if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] %s", timestamp, level, msg)
  cat(line, "\n", file = log_file, append = TRUE)
}

get_organism_db <- function(organism) {
  organism_lower <- tolower(organism)
  
  if (organism_lower %in% c("mouse", "mm", "mus musculus", "mus_musculus", "mmu")) {
    log_msg("Detected species: Mouse (Mus musculus)")
    return(list(
      CellChatDB = CellChatDB.mouse,
      PPI = PPI.mouse,
      species_name = "mouse"
    ))
  } else if (organism_lower %in% c("human", "hs", "homo sapiens", "homo_sapiens", "hsa")) {
    log_msg("Detected species: Human (Homo sapiens)")
    return(list(
      CellChatDB = CellChatDB.human,
      PPI = PPI.human,
      species_name = "human"
    ))
  } else {
    stop(paste("Unsupported species:", organism, "。Please use 'mouse' or 'human'"))
  }
}

is_seurat_v5 <- function() {
  seurat_version <- packageVersion("Seurat")
  return(seurat_version >= "5.0.0")
}

join_layers_if_needed <- function(seurat_obj, assay = "RNA") {
  if (!is_seurat_v5()) {
    log_msg("Detected Seurat V4, no need to join layers")
    return(seurat_obj)
  }
  
  log_msg("Seurat V5 detected, checking layers structure...")
  
  if (!assay %in% Assays(seurat_obj)) {
    log_msg(paste("warning：Assay", assay, "Does not exist"), level = "WARN")
    return(seurat_obj)
  }
  
  assay_obj <- seurat_obj[[assay]]
  assay_class <- class(assay_obj)[1]
  log_msg(paste("Assay type:", assay_class))
  
  if (assay_class == "SCTAssay") {
    log_msg("SCTAssay does not need JoinLayers, skipping")
    return(seurat_obj)
  }
  
  available_layers <- Layers(assay_obj)
  log_msg(paste("Available layers:", paste(available_layers, collapse = ", ")))
  
  counts_layers <- available_layers[grepl("^counts", available_layers)]
  data_layers <- available_layers[grepl("^data", available_layers)]
  
  need_join <- (length(counts_layers) > 1) || 
               (length(data_layers) > 1) ||
               (length(data_layers) == 1 && data_layers[1] != "data")
  
  if (need_join) {
    log_msg(paste("Detected dispersed layers, merging", assay, "assay..."))
    
    tryCatch({
      seurat_obj[[assay]] <- JoinLayers(seurat_obj[[assay]])
      log_msg(paste(assay, "assay layers merged successfully"))
    }, error = function(e) {
      log_msg(paste("JoinLayers failed:", e$message), level = "WARN")
    })
  }
  
  return(seurat_obj)
}

get_normalized_data <- function(seurat_obj, assay = "RNA") {
  if (is_seurat_v5()) {
    log_msg("Using Seurat V5 way to get expression matrix...")
    
    if (!is.null(attr(seurat_obj, "merged_data"))) {
      log_msg("Using manually merged data matrix")
      return(attr(seurat_obj, "merged_data"))
    }
    
    assay_obj <- seurat_obj[[assay]]
    available_layers <- Layers(assay_obj)
    
    if ("data" %in% available_layers) {
      expr <- LayerData(seurat_obj, assay = assay, layer = "data")
      log_msg("Get normalized data from the 'data' layer")
    } else {
      data_layers <- available_layers[grepl("^data", available_layers)]
      
      if (length(data_layers) > 0) {
        log_msg(paste("Merging dispersed 'data' layers:", paste(data_layers, collapse = ", ")))
        
        expr <- LayerData(seurat_obj, assay = assay, layer = data_layers[1])
        if (length(data_layers) > 1) {
          for (i in 2:length(data_layers)) {
            layer_data <- LayerData(seurat_obj, assay = assay, layer = data_layers[i])
            expr <- cbind(expr, layer_data)
          }
        }
        
        cell_order <- colnames(seurat_obj)
        if (all(cell_order %in% colnames(expr))) {
          expr <- expr[, cell_order]
        }
      } else if ("counts" %in% available_layers) {
        expr <- LayerData(seurat_obj, assay = assay, layer = "counts")
        log_msg("warning：Did not find 'data' layer, using 'counts' layer", level = "WARN")
      } else {
        first_layer <- available_layers[1]
        expr <- LayerData(seurat_obj, assay = assay, layer = first_layer)
        log_msg(paste("Using the first available layer:", first_layer))
      }
    }
  } else {
    log_msg("Using Seurat V4 way to get expression matrix...")
    expr <- GetAssayData(seurat_obj, assay = assay, slot = "data")
  }
  
  return(expr)
}

run_cellchat_single <- function(seurat_obj, group_name, organism_db, 
                                 cell_type_column, db_type, min_cells) {
  log_msg(paste("=== Processing group:", group_name, "==="))
  log_msg(paste("  Number of cells:", ncol(seurat_obj)))
  
  assay_to_use <- if ("RNA" %in% Assays(seurat_obj)) "RNA" else if ("SCT" %in% Assays(seurat_obj)) "SCT" else Assays(seurat_obj)[1]
  DefaultAssay(seurat_obj) <- assay_to_use
  
  seurat_obj <- join_layers_if_needed(seurat_obj, assay = assay_to_use)
  
  data.input <- get_normalized_data(seurat_obj, assay = assay_to_use)
  log_msg(paste("  Using Assay:", assay_to_use))
  log_msg(paste("  Expression matrix dimensions:", nrow(data.input), "genes x", ncol(data.input), "cells"))
  
  if (!inherits(data.input, "dgCMatrix")) {
    data.input <- as(data.input, "dgCMatrix")
  }
  
  labels <- as.character(seurat_obj@meta.data[[cell_type_column]])
  labels <- gsub("[^[:alnum:]_.-]", "_", labels)
  names(labels) <- colnames(seurat_obj)
  
  meta <- data.frame(
    labels = labels,
    row.names = names(labels),
    stringsAsFactors = FALSE
  )
  
  cell_types <- unique(labels)
  log_msg(paste("  Number of cell types:", length(cell_types)))
  
  log_msg("  Creating CellChat object...")
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  
  CellChatDB <- organism_db$CellChatDB
  
  if (!is.null(db_type) && db_type != "all" && db_type != "") {
    log_msg(paste("  Using database subset:", db_type))
    CellChatDB.use <- subsetDB(CellChatDB, search = db_type)
  } else {
    CellChatDB.use <- CellChatDB
  }
  
  cellchat@DB <- CellChatDB.use
  log_msg(paste("  Database contains", nrow(CellChatDB.use$interaction), "ligand-receptor pairs"))
  
  log_msg("  Subsetting data...")
  cellchat <- subsetData(cellchat)
  
  log_msg("  Identify overexpressed genes...")
  cellchat <- identifyOverExpressedGenes(cellchat)
  
  log_msg("  Identify overexpressed ligand-receptor interactions...")
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  log_msg("  Projecting to protein-protein interaction network...")
  cellchat <- projectData(cellchat, organism_db$PPI)
  
  log_msg("  Computing intercellular communication probabilities...")
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  
  log_msg(paste("  Filtering communication (minimum cells:", min_cells, "):"))
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  
  log_msg("  Computing pathway-level communication...")
  cellchat <- computeCommunProbPathway(cellchat)
  
  log_msg("  Aggregating cell communication networks...")
  cellchat <- aggregateNet(cellchat)
  
  log_msg(paste("  Group", group_name, "processing complete"))
  
  return(cellchat)
}

tryCatch({
  
  log_msg("========== Cell communication analysis begins==========")
  log_msg(paste("input file:", basename(input_file)))
  log_msg(paste("Seurat version:", as.character(packageVersion("Seurat"))))
  log_msg(paste("CellChat version:", as.character(packageVersion("CellChat"))))
  
  organism_db <- get_organism_db(organism)
  log_msg(paste("species:", organism_db$species_name))
  log_msg(paste("database type:", db_type))
  log_msg(paste("minimum cells:", min_cells))
  log_msg(paste("cell type column:", cell_type_column))
  
  for (dir_path in c(dirname(output_cellchat), dirname(output_net_df))) {
    if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
  }
  
  log_msg("Loading Seurat object...")
  seurat_obj <- readRDS(input_file)
  log_msg(paste("Loading complete: Contains", ncol(seurat_obj), "cells,", nrow(seurat_obj), "genes"))
  
  if (!cell_type_column %in% colnames(seurat_obj@meta.data)) {
    available_cols <- colnames(seurat_obj@meta.data)
    log_msg(paste("Warning: Specified cell type column not found:", cell_type_column), level = "WARN")
    
    possible_cols <- c("celltype_fine", "singler_labels", "cell_type", "celltype", "CellType",
                       "annotation", "seurat_clusters", "cluster")
    found <- FALSE
    for (col in possible_cols) {
      if (col %in% available_cols) {
        cell_type_column <- col
        log_msg(paste("Auto-selecting cell type column:", cell_type_column))
        found <- TRUE
        break
      }
    }
    if (!found) {
      stop("Unable to find a suitable cell type column")
    }
  }
  
  n_celltypes_overall <- length(unique(as.character(seurat_obj@meta.data[[cell_type_column]])))
  if (n_celltypes_overall < 2) {
    msg <- paste("CellChat skipped: Insufficient cell types (", n_celltypes_overall, ")", sep = "")
    log_msg(msg, level = "WARN")
    placeholder <- list(
      skipped = TRUE,
      reason = msg,
      n_celltypes = n_celltypes_overall
    )
    saveRDS(placeholder, output_cellchat)
    write.csv(data.frame(message = msg), output_net_df, row.names = FALSE)
    write.csv(data.frame(message = msg), output_lr_pairs, row.names = FALSE)
    log_msg("The placeholder CellChat results have been output.")
    quit(save = "no", status = 0)
  }
  
  has_groups <- "group" %in% colnames(seurat_obj@meta.data)
  
  if (has_groups) {
    groups <- unique(seurat_obj@meta.data$group)
    n_groups <- length(groups)
    log_msg(paste("Detected", n_groups, "groups:", paste(groups, collapse = ", ")))
    
    if (n_groups >= 2) {
      log_msg("Will perform CellChat analysis by group...")
      
      cellchat_list <- list()
      all_net_df <- data.frame()
      all_lr_pairs <- data.frame()
      
      for (grp in groups) {
        log_msg(paste("\nProcessing group:", grp))
        
        seurat_subset <- subset(seurat_obj, group == grp)
        
        if (ncol(seurat_subset) < 50) {
          log_msg(paste("Warning: Group", grp, "has too few cells (", ncol(seurat_subset), "), skipping"), level = "WARN")
          next
        }
        
        cellchat <- run_cellchat_single(
          seurat_obj = seurat_subset,
          group_name = grp,
          organism_db = organism_db,
          cell_type_column = cell_type_column,
          db_type = db_type,
          min_cells = min_cells
        )
        
        cellchat_list[[grp]] <- cellchat
        
        df.net <- subsetCommunication(cellchat)
        if (nrow(df.net) > 0) {
          df.net$group <- grp
          all_net_df <- rbind(all_net_df, df.net)
        }
        
        pathways <- cellchat@netP$pathways
        if (length(pathways) > 0) {
          for (pathway in pathways) {
            tryCatch({
              lr <- extractEnrichedLR(cellchat, signaling = pathway, geneLR.return = TRUE)
              if (nrow(lr) > 0) {
                lr$pathway <- pathway
                lr$group <- grp
                all_lr_pairs <- rbind(all_lr_pairs, lr)
              }
            }, error = function(e) {})
          }
        }
        
        group_output_dir <- file.path(dirname(output_cellchat), "by_group")
        if (!dir.exists(group_output_dir)) dir.create(group_output_dir, recursive = TRUE)
        saveRDS(cellchat, file.path(group_output_dir, paste0("cellchat_", grp, ".rds")))
        log_msg(paste("CellChat object for group", grp, "has been saved"))
      }
      
      if (length(cellchat_list) >= 2) {
        log_msg("\nMerging CellChat objects for comparative analysis...")
        
        cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
        log_msg("CellChat objects merged successfully")
        
        saveRDS(cellchat_merged, output_cellchat)
        log_msg(paste("Merged CellChat object saved to:", basename(output_cellchat)))
        
        saveRDS(cellchat_list, file.path(dirname(output_cellchat), "cellchat_list.rds"))
        
      } else if (length(cellchat_list) == 1) {
        log_msg("There is only one valid group, saving a single CellChat object")
        saveRDS(cellchat_list[[1]], output_cellchat)
      } else {
        stop("There are not enough valid groups for analysis")
      }
      
      if (nrow(all_net_df) > 0) {
        write.csv(all_net_df, output_net_df, row.names = FALSE, quote = FALSE)
        log_msg(paste("Communication network has been saved (", nrow(all_net_df), " records)"))
      } else {
        write.csv(data.frame(message = "No significant interactions detected"), 
                  output_net_df, row.names = FALSE)
      }
      
      if (nrow(all_lr_pairs) > 0) {
        write.csv(all_lr_pairs, output_lr_pairs, row.names = FALSE, quote = FALSE)
        log_msg(paste("Ligand-receptor pair information has been saved (", nrow(all_lr_pairs), " records)"))
      } else {
        write.csv(data.frame(message = "No ligand-receptor pairs extracted"), 
                  output_lr_pairs, row.names = FALSE)
      }
      
    } else {
      log_msg("There is only one group, performing overall analysis...")
      has_groups <- FALSE
    }
  }
  
  if (!has_groups) {
    log_msg("Performing overall CellChat analysis (no groups)...")
    
    cellchat <- run_cellchat_single(
      seurat_obj = seurat_obj,
      group_name = "all",
      organism_db = organism_db,
      cell_type_column = cell_type_column,
      db_type = db_type,
      min_cells = min_cells
    )
    
    df.net <- subsetCommunication(cellchat)
    
    if (nrow(df.net) > 0) {
      log_msg(paste("Significant ligand-receptor interactions detected:", nrow(df.net)))
      write.csv(df.net, output_net_df, row.names = FALSE, quote = FALSE)
    } else {
      log_msg("Warning: No significant ligand-receptor interactions detected", level = "WARN")
      write.csv(data.frame(message = "No significant interactions detected"), 
                output_net_df, row.names = FALSE)
    }
    
    pathways <- cellchat@netP$pathways
    
    if (!is.null(pathways) && length(pathways) > 0) {
      log_msg(paste("Signaling pathways detected:", length(pathways)))
      
      all_lr_pairs <- data.frame()
      for (pathway in pathways) {
        tryCatch({
          lr <- extractEnrichedLR(cellchat, signaling = pathway, geneLR.return = TRUE)
          if (nrow(lr) > 0) {
            lr$pathway <- pathway
            all_lr_pairs <- rbind(all_lr_pairs, lr)
          }
        }, error = function(e) {})
      }
      
      if (nrow(all_lr_pairs) > 0) {
        write.csv(all_lr_pairs, output_lr_pairs, row.names = FALSE, quote = FALSE)
        log_msg(paste("Ligand-receptor pair information has been saved (", nrow(all_lr_pairs), " records)"))
      } else {
        write.csv(data.frame(message = "No ligand-receptor pairs extracted"), 
                  output_lr_pairs, row.names = FALSE)
      }
    } else {
      log_msg("Warning: No signaling pathways detected", level = "WARN")
      write.csv(data.frame(message = "No signaling pathways detected"), 
                output_lr_pairs, row.names = FALSE)
    }
    
    saveRDS(cellchat, output_cellchat)
    log_msg(paste("CellChat object has been saved to:", basename(output_cellchat)))
  }
  
  summary_text <- paste(
    "",
    "==========================================",
    "  Summary of CellChat Inference Analysis",
    "==========================================",
    paste("Species:", organism_db$species_name),
    paste("Group Analysis:", ifelse(has_groups && n_groups >= 2, paste("Yes (", n_groups, " groups)"), "No")),
    paste("Database Used:", db_type),
    paste("Assay Used: RNA (log-normalized)"),
    "==========================================",
    "",
    sep = "\n"
  )
  cat(summary_text)
  log_msg("CellChat Communication inference analysis completed!")
  
}, error = function(e) {
  log_msg(paste("error:", e$message), level = "ERROR")
  log_msg(paste("Error traceback:", paste(capture.output(traceback()), collapse = "\n")), level = "ERROR")
  msg <- paste("CellChat inference failed, outputting placeholder results:", e$message)
  placeholder <- list(skipped = TRUE, reason = msg)
  saveRDS(placeholder, output_cellchat)
  write.csv(data.frame(message = msg), output_net_df, row.names = FALSE)
  write.csv(data.frame(message = msg), output_lr_pairs, row.names = FALSE)
  log_msg(msg, level = "WARN")
  quit(save = "no", status = 0)
})
