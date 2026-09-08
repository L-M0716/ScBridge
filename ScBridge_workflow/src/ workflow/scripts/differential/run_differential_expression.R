# scripts/differential_scripts/run_differential_expression.R

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(DESeq2)
  library(Matrix)
  library(tools)
})

input_seurat       <- snakemake@input[["annotated_seurat"]]
output_csv         <- snakemake@output[["de_results_table"]]
log_file           <- snakemake@log[[1]]
threads            <- snakemake@threads
samples_table_path <- snakemake@input[["sample_table"]]

group_by_col    <- snakemake@params[["group_by_column"]]
control_group   <- snakemake@params[["control_group"]]
logfc_threshold <- snakemake@params[["logfc_threshold"]]
min_cells_group <- snakemake@params[["min_cells_group"]]

if (!dir.exists(dirname(log_file))) dir.create(dirname(log_file), recursive = TRUE)
log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] %s", timestamp, level, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  message(line)
}

tryCatch({
  log_msg("Differential expression analysis.")

  ext <- tolower(file_ext(samples_table_path))
  sample_info <- if (ext == "csv") {
    read.csv(samples_table_path, stringsAsFactors = FALSE)
  } else {
    read.delim(samples_table_path, stringsAsFactors = FALSE)
  }
  log_msg(paste("samples:", nrow(sample_info), "samples"))
  log_msg(paste("sample IDs:", paste(sample_info$sample_id, collapse = ", ")))
  log_msg(paste("groups:", paste(sample_info[[group_by_col]], collapse = ", ")))

  if (!group_by_col %in% colnames(sample_info))
    stop(paste("Grouping column not found in the sample table:", group_by_col))

  all_groups <- unique(sample_info[[group_by_col]])
  if (!control_group %in% all_groups)
    stop(paste("Sample table does not contain the control group:", control_group))

  treated_group <- if (length(all_groups) == 2) {
    setdiff(all_groups, control_group)
  } else {
    g <- setdiff(all_groups, control_group)[1]
    log_msg(paste("Multiple groups found, using the first non-control group:", g), level = "WARN")
    g
  }
  log_msg(paste("Compare:", treated_group, "(Handle) vs", control_group, "(对照)"))

  log_msg(paste("loading Seurat object:", basename(input_seurat)))
  seurat_obj <- readRDS(input_seurat)
  available_assays <- names(seurat_obj@assays)
  log_msg(paste("available assays:", paste(available_assays, collapse = ", ")))

  if ("RNA" %in% available_assays) {
    use_assay <- "RNA"
    log_msg("Using RNA assay (log-normalized counts, recommended for DE）")
  } else if ("SCT" %in% available_assays) {
    use_assay <- "SCT"
    log_msg("Using SCT assay", level = "WARN")
  } else {
    stop("No assay available for differential analysis (RNA or SCT).")
  }
  DefaultAssay(seurat_obj) <- use_assay

  if (!group_by_col %in% colnames(seurat_obj@meta.data)) {
    log_msg("Seurat object does not contain a group column, adding from the sample table...")
    sample_col_meta <- if ("sample_id" %in% colnames(seurat_obj@meta.data)) {
      "sample_id"
    } else {
      "orig.ident"
    }
    sample_to_group <- setNames(sample_info[[group_by_col]], sample_info$sample_id)
    seurat_obj@meta.data[[group_by_col]] <- sample_to_group[seurat_obj@meta.data[[sample_col_meta]]]
    na_n <- sum(is.na(seurat_obj@meta.data[[group_by_col]]))
    if (na_n > 0) {
      log_msg(paste("warning:", na_n, "Cells that cannot be matched to a group will be removed"), level = "WARN")
      seurat_obj <- subset(seurat_obj,
        cells = rownames(seurat_obj@meta.data)[!is.na(seurat_obj@meta.data[[group_by_col]])])
    }
    log_msg("group column added.")
  }

  group_table <- table(seurat_obj@meta.data[[group_by_col]])
  for (g in names(group_table))
    log_msg(paste("  ", g, ":", group_table[g], "cells"))

  celltype_col <- if ("celltype_fine" %in% colnames(seurat_obj@meta.data)) {
    "celltype_fine"
  } else if ("singler_labels" %in% colnames(seurat_obj@meta.data)) {
    "singler_labels"
  } else if ("manual_annotation" %in% colnames(seurat_obj@meta.data)) {
    "manual_annotation"
  } else {
    stop("No cell type annotation column found (celltype_fine / singler_labels / manual_annotation).")
  }
  log_msg(paste("Using cell type annotation column:", celltype_col))

  sample_col <- if ("sample_id" %in% colnames(seurat_obj@meta.data)) "sample_id" else "orig.ident"

  all_cell_types <- sort(unique(na.omit(seurat_obj@meta.data[[celltype_col]])))
  all_de_results <- list()
  log_msg(paste("Cell types:", paste(all_cell_types, collapse = ", ")))

  for (cell_type in all_cell_types) {
    log_msg(paste("\n--- Cell type:", cell_type, "---"))

    cells_of_type <- rownames(seurat_obj@meta.data)[
      seurat_obj@meta.data[[celltype_col]] == cell_type]
    if (length(cells_of_type) == 0) {
      log_msg("  No cells found, skipping", level = "WARN"); next
    }

    seurat_sub <- subset(seurat_obj, cells = cells_of_type)
    meta_sub   <- seurat_sub@meta.data

    group_counts <- table(meta_sub[[group_by_col]])
    log_msg(paste("  Cell distribution:", paste(names(group_counts), "=", group_counts, collapse = ", ")))

    if (!control_group %in% names(group_counts) || !treated_group %in% names(group_counts)) {
      log_msg("  Required groups missing, skipping", level = "WARN"); next
    }
    if (group_counts[control_group] < min_cells_group ||
        group_counts[treated_group] < min_cells_group) {
      log_msg(paste("  Not enough cells (min =", min_cells_group, "), skipping"), level = "WARN"); next
    }

    samples_here <- unique(meta_sub[[sample_col]])
    counts_mat   <- GetAssayData(seurat_sub, assay = use_assay, layer = "counts")
    if (is.null(counts_mat) || nrow(counts_mat) == 0) {
      log_msg("  counts layer is empty, skipping", level = "WARN"); next
    }

    pb_list <- lapply(samples_here, function(sid) {
      idx <- which(meta_sub[[sample_col]] == sid)
      if (length(idx) == 0) return(NULL)
      Matrix::rowSums(counts_mat[, idx, drop = FALSE])
    })
    names(pb_list) <- samples_here
    pb_list <- pb_list[!sapply(pb_list, is.null)]

    pb_mat <- do.call(cbind, pb_list)

    sample_group_map <- unique(meta_sub[, c(sample_col, group_by_col)])
    col_data <- data.frame(
      sample_id = colnames(pb_mat),
      row.names  = colnames(pb_mat),
      stringsAsFactors = FALSE
    )
    col_data[[group_by_col]] <- sample_group_map[[group_by_col]][
      match(col_data$sample_id, sample_group_map[[sample_col]])]
    col_data[[group_by_col]] <- factor(col_data[[group_by_col]],
                                       levels = c(control_group, treated_group))
    col_data <- col_data[!is.na(col_data[[group_by_col]]), , drop = FALSE]
    pb_mat   <- pb_mat[, rownames(col_data), drop = FALSE]

    n_ctrl  <- sum(col_data[[group_by_col]] == control_group)
    n_treat <- sum(col_data[[group_by_col]] == treated_group)
    log_msg(paste("  Pseudo-bulk samples: ctrl =", n_ctrl, ", treat =", n_treat))

    if (n_ctrl < 2 || n_treat < 2) {
      log_msg("  Each group needs at least 2 samples to run DESeq2, skipping", level = "WARN"); next
    }

    keep   <- rowSums(pb_mat >= 1) >= min(n_ctrl, n_treat)
    pb_mat <- pb_mat[keep, , drop = FALSE]
    log_msg(paste("  Genes after filtering:", nrow(pb_mat)))

    de_res <- tryCatch({
      dds <- DESeqDataSetFromMatrix(
        countData = round(as.matrix(pb_mat)),
        colData   = col_data,
        design    = as.formula(paste("~", group_by_col))
      )
      dds <- DESeq(dds, quiet = TRUE)
      res <- results(dds,
                     contrast     = c(group_by_col, treated_group, control_group),
                     alpha        = 0.05,
                     lfcThreshold = logfc_threshold)
      res_df <- as.data.frame(res)
      res_df$gene       <- rownames(res_df)
      res_df$cell_type  <- cell_type
      res_df$comparison <- paste0(treated_group, "_vs_", control_group)
      res_df$assay_used <- use_assay
      res_df$avg_log2FC <- res_df$log2FoldChange
      res_df$p_val      <- res_df$pvalue
      res_df$p_val_adj  <- res_df$padj
      res_df[!is.na(res_df$p_val_adj), ]
    }, error = function(e) {
      log_msg(paste("  DESeq2 失败:", e$message), level = "ERROR")
      NULL
    })

    if (!is.null(de_res) && nrow(de_res) > 0) {
      sig   <- sum(de_res$p_val_adj < 0.05, na.rm = TRUE)
      up    <- sum(de_res$p_val_adj < 0.05 & de_res$avg_log2FC > 0, na.rm = TRUE)
      down  <- sum(de_res$p_val_adj < 0.05 & de_res$avg_log2FC < 0, na.rm = TRUE)
      log_msg(paste("  Genes detected:", nrow(de_res), "| Significant:", sig,
                    "( Upregulated:", up, ", Downregulated:", down, ")"))
      all_de_results[[cell_type]] <- de_res
    } else {
      log_msg("  No differentially expressed genes found")
    }
  }

  if (length(all_de_results) > 0) {
    final_df <- bind_rows(all_de_results)
    final_df$significant <- final_df$p_val_adj < 0.05
    log_msg(paste("\nCompletion. Total of", length(all_de_results), "cell types,",
                  nrow(final_df), "records."))
  } else {
    log_msg("  No differentially expressed genes found in any cell type.", level = "WARN")
    final_df <- data.frame(
      gene = character(), cell_type = character(), comparison = character(),
      avg_log2FC = numeric(), p_val = numeric(), p_val_adj = numeric(),
      significant = logical(), assay_used = character()
    )
  }

  if (!dir.exists(dirname(output_csv))) dir.create(dirname(output_csv), recursive = TRUE)
  write.csv(final_df, file = output_csv, row.names = FALSE)
  log_msg(paste("result:", basename(output_csv)))
  log_msg("Differential expression analysis workflow successfully completed！")

}, error = function(e) {
  log_msg(paste("FATAL ERROR:", e$message), level = "ERROR")
  stop(e)
})
