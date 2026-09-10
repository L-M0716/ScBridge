# scripts/monocle3/monocle3_analysis.R
Sys.setenv(ICU_CACHE_DATA = "TRUE")
suppressMessages({
  library(stringi)
  invisible(stringi::stri_info())
})

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(tidyverse)
  library(patchwork)
  library(Matrix)
  library(pheatmap)
  library(scales)
  library(viridis)
})
input_file       <- snakemake@input[["annotated"]]
output_cds       <- snakemake@output[["cds_object"]]
output_trajectory <- snakemake@output[["trajectory_summary"]]
output_genes     <- snakemake@output[["trajectory_genes"]]
output_plots_dir <- snakemake@output[["plot_dir"]]
clusters_to_use   <- snakemake@params[["clusters_to_use"]]
root_cell_type    <- snakemake@params[["root_cell_type"]]
root_cluster      <- snakemake@params[["root_cluster"]]
num_dim           <- snakemake@params[["num_dim"]]
resolution        <- snakemake@params[["resolution"]]
k                 <- snakemake@params[["k"]]
q_value_threshold <- snakemake@params[["q_value_threshold"]]
log_file <- snakemake@log[[1]]


.get_expr <- function(obj, assay, want_counts = TRUE) {
  layer_name <- if (want_counts) "counts" else "data"
  if (packageVersion("SeuratObject") >= "5.0.0") {
    Seurat::GetAssayData(obj, assay = assay, layer = layer_name)
  } else {
    Seurat::GetAssayData(obj, assay = assay, slot = layer_name)
  }
}

seurat_to_monocle3 <- function(seurat_obj, assay = "RNA", use_raw_counts = TRUE) {
  expr <- .get_expr(seurat_obj, assay, want_counts = use_raw_counts)
  rn <- rownames(expr)
  if (anyDuplicated(rn)) {
    warning("Duplicate gene names detected and automatically deduplicated (make.unique).")
    rownames(expr) <- make.unique(rn)
  }
  cell_md <- seurat_obj@meta.data
  gene_md <- data.frame(
    gene_short_name = rownames(expr),
    row.names = rownames(expr),
    stringsAsFactors = FALSE
  )
  
  cds <- monocle3::new_cell_data_set(
    expression_data = expr,
    cell_metadata = cell_md,
    gene_metadata = gene_md
  )
  
  if (!"seurat_clusters" %in% colnames(SummarizedExperiment::colData(cds))) {
    SummarizedExperiment::colData(cds)$seurat_clusters <- as.character(Seurat::Idents(seurat_obj))
  }
  cds
}

get_valid_root_pr_node <- function(cds) {
  g <- monocle3::principal_graph(cds)[["UMAP"]]
  if (is.null(g)) return(NA_character_)
  vnames <- igraph::V(g)$name
  
  aux <- monocle3::principal_graph_aux(cds)[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  
  nn <- NULL
  if (is.matrix(aux)) {
    nn <- aux[, 1]
  } else if (is.data.frame(aux)) {
    nn <- aux[[1]]
  } else if (is.list(aux) && length(aux) == 1) {
    x <- aux[[1]]
    if (is.matrix(x)) nn <- x[, 1]
    else if (is.data.frame(x)) nn <- x[[1]]
    else nn <- x
  } else {
    nn <- aux
  }
  nn <- as.character(nn)
  
  nn <- nn[nn %in% vnames]
  if (length(nn) == 0) return(NA_character_)
  names(sort(table(nn), decreasing = TRUE))[1]
}

if (!dir.exists(dirname(log_file))) dir.create(dirname(log_file), recursive = TRUE)
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] %s", timestamp, level, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  message(line)
}

custom_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "right",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )

save_png <- function(plot_obj, filename, width = 3000, height = 2500, res = 300) {
  png(filename, width = width, height = height, res = res)
  print(plot_obj)
  dev.off()
}

tryCatch({
  
  log_msg("Start Monocle3 pseudotime analysis workflow")
  log_msg(paste("Input file:", basename(input_file)))
  
  log_msg("Loading Seurat object...")
  seurat_obj <- readRDS(input_file)
  log_msg(paste("Loading completed, containing", ncol(seurat_obj), "cells, ", nrow(seurat_obj), "genes"))
  
  if (!is.null(clusters_to_use) && clusters_to_use != "all" && clusters_to_use != "") {
    clusters_list <- unlist(strsplit(clusters_to_use, ","))
    log_msg(paste("Filter the specified cluster:", paste(clusters_list, collapse = ", ")))
    seurat_obj <- subset(seurat_obj, subset = seurat_clusters %in% clusters_list)
    log_msg(paste("Filtered, remaining", ncol(seurat_obj), "cells"))
  }
  
  assay_names <- names(seurat_obj@assays)
  assay_for_monocle <- if ("RNA" %in% assay_names) "RNA" else if ("SCT" %in% assay_names) "SCT" else assay_names[1]
  log_msg(paste("Convert Seurat object to CDS object (without using SeuratWrappers), assay =", assay_for_monocle))
  cds <- seurat_to_monocle3(seurat_obj, assay = assay_for_monocle, use_raw_counts = TRUE)
  
  log_msg("Estimate cell size factors...")
  cds <- estimate_size_factors(cds)
  
  if (!"gene_short_name" %in% colnames(SummarizedExperiment::rowData(cds))) {
    SummarizedExperiment::rowData(cds)$gene_short_name <- rownames(cds)
  }
  SummarizedExperiment::rowData(cds)$gene_name <- rownames(cds)
  
  celltype_source <- if ("celltype_fine" %in% colnames(seurat_obj@meta.data)) "celltype_fine" else
                     if ("singler_labels" %in% colnames(seurat_obj@meta.data)) "singler_labels" else
                     if ("manual_annotation" %in% colnames(seurat_obj@meta.data)) "manual_annotation" else NULL
  if (!is.null(celltype_source)) {
    SummarizedExperiment::colData(cds)$cell_type <- seurat_obj@meta.data[[celltype_source]]
    log_msg(paste("Transferred cell type annotation information, source column:", celltype_source))
  }
  if (!"seurat_clusters" %in% colnames(SummarizedExperiment::colData(cds))) {
    SummarizedExperiment::colData(cds)$seurat_clusters <- seurat_obj$seurat_clusters
  }
  if ("sample_id" %in% colnames(seurat_obj@meta.data) &&
      !"sample_id" %in% colnames(SummarizedExperiment::colData(cds))) {
    SummarizedExperiment::colData(cds)$sample_id <- seurat_obj$sample_id
  }
  if ("group" %in% colnames(seurat_obj@meta.data) &&
      !"group" %in% colnames(SummarizedExperiment::colData(cds))) {
    SummarizedExperiment::colData(cds)$group <- seurat_obj$group
  }
  
  log_msg(paste("Executing CDS preprocessing, dimensions:", num_dim))
  cds <- preprocess_cds(
    cds,
    num_dim = num_dim,
    norm_method = "size_only",
    verbose = TRUE
  )
  
  log_msg("Executing UMAP dimensionality reduction...")
  cds <- reduce_dimension(
    cds,
    max_components = 2,
    reduction_method = "UMAP",
    preprocess_method = "PCA",
    verbose = TRUE
  )
  
  log_msg(paste("Performing cell clustering, resolution:", resolution, ", k-value:", k))
  cds <- cluster_cells(
    cds,
    resolution = resolution,
    k = k,
    verbose = TRUE
  )
  
  log_msg("Learning cell trajectory...")
  cds <- learn_graph(cds, verbose = TRUE)
  
  log_msg("Ordering cells in pseudotime...")
  
  get_root_cells <- function(cds, root_cell_type = NULL, root_cluster = NULL) {
    root_cells <- NULL
    if (!is.null(root_cell_type) && root_cell_type != "" &&
        "cell_type" %in% colnames(SummarizedExperiment::colData(cds))) {
      root_cells <- colnames(cds)[SummarizedExperiment::colData(cds)$cell_type == root_cell_type]
      log_msg(paste("Using cell type", root_cell_type, "as root node, containing", length(root_cells), "cells"))
    } else if (!is.null(root_cluster) && root_cluster != "") {
      root_cells <- colnames(cds)[SummarizedExperiment::colData(cds)$seurat_clusters == root_cluster]
      log_msg(paste("Using cluster", root_cluster, "as root node, containing", length(root_cells), "cells"))
    }
    root_cells
  }
  
  root_cells <- get_root_cells(cds, root_cell_type, root_cluster)
  
  if (!is.null(root_cells) && length(root_cells) > 0) {
    cds <- order_cells(cds, root_cells = root_cells)
    log_msg("Using manually specified root cells to complete pseudotime ordering")
  } else {
    log_msg("Automatically selecting root node...")
    root_node <- get_valid_root_pr_node(cds)
    g <- monocle3::principal_graph(cds)[["UMAP"]]
    vnames <- if (!is.null(g)) igraph::V(g)$name else character(0)
    
    if (!is.na(root_node) && root_node %in% vnames) {
      cds <- order_cells(cds, root_pr_nodes = root_node)
      log_msg(paste("Using automatically selected root node to complete pseudotime ordering:", root_node))
    } else {
      log_msg("Failed to reliably identify root_pr_node, switching to root cell strategy (based on UMAP1's starting cluster)")
      umap <- SingleCellExperiment::reducedDims(cds)$UMAP
      if (is.null(umap)) stop("UMAP embedding not found, cannot automatically select root cells")
      clu <- as.character(SummarizedExperiment::colData(cds)$seurat_clusters)
      if (is.null(clu)) stop("seurat_clusters not found, cannot automatically select root cells")
      
      sc <- tapply(umap[, 1], clu, median, na.rm = TRUE)
      root_cluster_auto <- names(which.min(sc))
      root_cells_auto <- colnames(cds)[clu == root_cluster_auto]
      root_cells_auto <- head(root_cells_auto, 500)  # 可调整
      
      cds <- order_cells(cds, root_cells = root_cells_auto)
      log_msg(paste("Using cluster", root_cluster_auto, "with", length(root_cells_auto), "cells as root cells for pseudotime ordering"))
    }
  }
  
  SummarizedExperiment::colData(cds)$monocle3_pseudotime <- monocle3::pseudotime(cds)
  log_msg("Pseudotime calculation complete")
  
  log_msg("Creating visualization results...")
  if (!dir.exists(output_plots_dir)) {
    dir.create(output_plots_dir, recursive = TRUE)
  }
  
  n_clusters <- length(unique(SummarizedExperiment::colData(cds)$seurat_clusters))
  has_celltype <- "cell_type" %in% colnames(SummarizedExperiment::colData(cds))
  n_celltypes <- if(has_celltype) length(unique(SummarizedExperiment::colData(cds)$cell_type)) else 0
  has_group <- "group" %in% colnames(SummarizedExperiment::colData(cds))
  
  log_msg(paste("Detected", n_clusters, "clusters,", n_celltypes, "cell types"))
  
  pdf(file.path(output_plots_dir, "trajectory_main.pdf"), 
      width = 14, height = 12)
  
  p1 <- plot_cells(cds,
                   color_cells_by = "seurat_clusters",
                   label_cell_groups = TRUE,
                   label_leaves = FALSE,
                   label_branch_points = FALSE,
                   group_label_size = 4,
                   cell_size = 0.5,
                   alpha = 0.8) +
    ggtitle("Trajectory Colored by Clusters") +
    custom_theme +
    guides(color = guide_legend(ncol = ifelse(n_clusters > 15, 2, 1),
                                override.aes = list(size = 3, alpha = 1)))
  print(p1)
  
  p1b <- plot_cells(cds,
                    color_cells_by = "seurat_clusters",
                    label_cell_groups = FALSE,
                    label_leaves = TRUE,
                    label_branch_points = TRUE,
                    graph_label_size = 3,
                    cell_size = 0.3,
                    alpha = 0.6) +
    ggtitle("Trajectory Structure",
            subtitle = "Showing branch points (black) and leaf nodes") +
    custom_theme +
    guides(color = guide_legend(ncol = ifelse(n_clusters > 15, 2, 1),
                                override.aes = list(size = 3, alpha = 1)))
  print(p1b)
  
  if (has_celltype) {
    p2 <- plot_cells(cds,
                     color_cells_by = "cell_type",
                     label_cell_groups = (n_celltypes <= 12),
                     label_leaves = FALSE,
                     label_branch_points = FALSE,
                     group_label_size = 3.5,
                     cell_size = 0.5,
                     alpha = 0.8) +
      ggtitle("Trajectory Colored by Cell Type") +
      custom_theme +
      guides(color = guide_legend(ncol = ifelse(n_celltypes > 10, 2, 1),
                                  override.aes = list(size = 3, alpha = 1)))
    print(p2)
  }
  
  p3 <- plot_cells(cds,
                   color_cells_by = "pseudotime",
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   label_branch_points = FALSE,
                   cell_size = 0.5,
                   alpha = 0.8) +
    ggtitle("Trajectory Colored by Pseudotime") +
    custom_theme +
    scale_color_viridis_c(option = "plasma", name = "Pseudotime")
  print(p3)
  
  p3b <- plot_cells(cds,
                    color_cells_by = "pseudotime",
                    label_cell_groups = FALSE,
                    label_leaves = TRUE,
                    label_branch_points = TRUE,
                    graph_label_size = 3,
                    cell_size = 0.4,
                    alpha = 0.7) +
    ggtitle("Pseudotime with Trajectory Graph",
            subtitle = "Branch points and leaf nodes labeled") +
    custom_theme +
    scale_color_viridis_c(option = "plasma", name = "Pseudotime")
  print(p3b)
  
  p4 <- plot_cells(cds,
                   color_cells_by = "partition",
                   label_cell_groups = TRUE,
                   show_trajectory_graph = FALSE,
                   group_label_size = 5,
                   cell_size = 0.5,
                   alpha = 0.8) +
    ggtitle("Cell Partitions") +
    custom_theme
  print(p4)
  
  if (has_group) {
    p5 <- plot_cells(cds,
                     color_cells_by = "group",
                     label_cell_groups = TRUE,
                     label_leaves = FALSE,
                     label_branch_points = FALSE,
                     group_label_size = 4,
                     cell_size = 0.5,
                     alpha = 0.8) +
      ggtitle("Trajectory Colored by Experimental Group") +
      custom_theme
    print(p5)
  }
  
  p_combined <- (p1 + p3) / (p1b + p4) +
    plot_annotation(
      title = "Monocle3 Trajectory Analysis Overview",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  print(p_combined)
  
  dev.off()
  log_msg(paste("Main trajectory plot saved to:", file.path(output_plots_dir, "trajectory_main.pdf")))
  
  log_msg("Generating high-resolution PNG images...")
  
  save_png(p1, file.path(output_plots_dir, "trajectory_clusters.png"))
  
  save_png(p3, file.path(output_plots_dir, "trajectory_pseudotime.png"))
  
  save_png(p1b, file.path(output_plots_dir, "trajectory_structure.png"))
  
  if (has_celltype) {
    save_png(p2, file.path(output_plots_dir, "trajectory_celltype.png"))
  }
  
  log_msg("High-resolution PNG images saved")
  
  log_msg("Analyzing pseudotime distribution...")
  
  pdf(file.path(output_plots_dir, "pseudotime_analysis.pdf"), width = 12, height = 10)
  
  data_pseudo <- as.data.frame(SummarizedExperiment::colData(cds))
  
  p_box <- ggplot(data_pseudo, aes(x = reorder(seurat_clusters, monocle3_pseudotime),
                                   y = monocle3_pseudotime,
                                   fill = seurat_clusters)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.4) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      legend.position = "none"
    ) +
    labs(x = "Cluster", y = "Pseudotime", 
         title = "Pseudotime Distribution by Cluster") +
    scale_fill_viridis_d(option = "turbo")
  print(p_box)
  
  p_violin <- ggplot(data_pseudo, aes(x = reorder(seurat_clusters, monocle3_pseudotime),
                                      y = monocle3_pseudotime,
                                      fill = seurat_clusters)) +
    geom_violin(scale = "width", alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.2) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      legend.position = "none"
    ) +
    labs(x = "Cluster", y = "Pseudotime", 
         title = "Pseudotime Distribution (Violin Plot)") +
    scale_fill_viridis_d(option = "turbo")
  print(p_violin)
  
  p_density <- ggplot(data_pseudo, aes(x = monocle3_pseudotime, fill = seurat_clusters)) +
    geom_density(alpha = 0.4) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      legend.text = element_text(size = 8),
      legend.position = "right"
    ) +
    labs(x = "Pseudotime", y = "Density", 
         title = "Pseudotime Density Distribution",
         fill = "Cluster") +
    scale_fill_viridis_d(option = "turbo") +
    guides(fill = guide_legend(ncol = ifelse(n_clusters > 12, 2, 1),
                               override.aes = list(alpha = 0.8)))
  print(p_density)
  
  if (n_clusters > 6) {
    p_density_facet <- ggplot(data_pseudo, aes(x = monocle3_pseudotime, fill = seurat_clusters)) +
      geom_density(alpha = 0.8) +
      geom_rug(alpha = 0.1, length = unit(0.02, "npc")) +
      facet_wrap(~seurat_clusters, scales = "free_y", ncol = 4) +
      theme_minimal() +
      theme(
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "gray95", color = NA),
        legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")
      ) +
      labs(x = "Pseudotime", y = "Density", 
           title = "Pseudotime Density by Cluster (Faceted)") +
      scale_fill_viridis_d(option = "turbo")
    print(p_density_facet)
  }
  
  if (has_celltype) {
    p_box_celltype <- ggplot(data_pseudo, aes(x = reorder(cell_type, monocle3_pseudotime),
                                              y = monocle3_pseudotime,
                                              fill = cell_type)) +
      geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.4) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.position = "none"
      ) +
      labs(x = "Cell Type", y = "Pseudotime", 
           title = "Pseudotime Distribution by Cell Type") +
      scale_fill_viridis_d(option = "mako")
    print(p_box_celltype)
    
    p_density_celltype <- ggplot(data_pseudo, aes(x = monocle3_pseudotime, fill = cell_type)) +
      geom_density(alpha = 0.5) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.position = "right"
      ) +
      labs(x = "Pseudotime", y = "Density", 
           title = "Pseudotime Density by Cell Type",
           fill = "Cell Type") +
      scale_fill_viridis_d(option = "mako") +
      guides(fill = guide_legend(ncol = ifelse(n_celltypes > 10, 2, 1)))
    print(p_density_celltype)
  }
  
  if (has_group) {
    p_box_group <- ggplot(data_pseudo, aes(x = group,
                                           y = monocle3_pseudotime,
                                           fill = group)) +
      geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
      geom_jitter(width = 0.2, alpha = 0.1, size = 0.3) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.position = "none"
      ) +
      labs(x = "Group", y = "Pseudotime", 
           title = "Pseudotime Distribution by Experimental Group")
    print(p_box_group)
    
    p_group_cluster <- ggplot(data_pseudo, aes(x = seurat_clusters,
                                               y = monocle3_pseudotime,
                                               fill = group)) +
      geom_boxplot(outlier.size = 0.2, position = position_dodge(0.8)) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold")
      ) +
      labs(x = "Cluster", y = "Pseudotime", 
           title = "Pseudotime by Cluster and Group",
           fill = "Group")
    print(p_group_cluster)
  }
  
  data_pseudo_sorted <- data_pseudo[order(data_pseudo$monocle3_pseudotime), ]
  p_scatter <- ggplot(data_pseudo_sorted, aes(x = seq_len(nrow(data_pseudo_sorted)),
                                       y = monocle3_pseudotime,
                                       color = monocle3_pseudotime)) +
    geom_point(size = 0.3, alpha = 0.5) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      legend.position = "right"
    ) +
    labs(x = "Cell Rank", y = "Pseudotime", 
         title = "Sorted Pseudotime Distribution") +
    scale_color_viridis_c(option = "plasma", name = "Pseudotime")
  print(p_scatter)
  
  dev.off()
  log_msg(paste("The planned timing analysis diagram has been saved to:", file.path(output_plots_dir, "pseudotime_analysis.pdf")))
  
  log_msg("Identify trajectory-dependent genes...")
  
  pr_test_res <- tryCatch({
    graph_test(cds,
               neighbor_graph = "principal_graph",
               cores = 1)
  }, error = function(e) {
    log_msg(paste("graph_test failed, falling back to empty result:", e$message), level = "WARN")
    data.frame()
  })
  
  if (nrow(pr_test_res) > 0 && all(c("q_value", "morans_I") %in% colnames(pr_test_res))) {
    pr_deg <- pr_test_res %>%
      dplyr::filter(q_value < q_value_threshold) %>%
      dplyr::arrange(dplyr::desc(morans_I))
  } else {
    pr_deg <- data.frame()
  }
  
  log_msg(paste("Found", nrow(pr_deg), "trajectory-dependent genes (q-value <", q_value_threshold, ")"))
  
  write.csv(pr_deg, output_genes, row.names = FALSE)
  log_msg(paste("Trajectory-dependent genes saved to:", basename(output_genes)))
  
  if (nrow(pr_deg) > 0) {
    log_msg("Visualizing top trajectory-dependent genes...")
    
    pdf(file.path(output_plots_dir, "gene_dynamics.pdf"), width = 16, height = 14)
    
    top_genes_9 <- head(rownames(pr_deg), 9)
    log_msg(paste("Draw Top", length(top_genes_9), "trajectory-dependent genes"))
    
    p_genes_umap <- plot_cells(cds,
                               genes = top_genes_9,
                               show_trajectory_graph = TRUE,
                               label_cell_groups = FALSE,
                               label_leaves = FALSE,
                               label_branch_points = FALSE,
                               cell_size = 0.15,
                               alpha = 0.6) +
      theme(
        strip.text = element_text(size = 11, face = "bold"),
        strip.background = element_rect(fill = "gray95", color = NA)
      )
    print(p_genes_umap)
    
    p_genes_clean <- plot_cells(cds,
                                genes = top_genes_9,
                                show_trajectory_graph = FALSE,
                                label_cell_groups = FALSE,
                                cell_size = 0.2,
                                alpha = 0.7) +
      theme(
        strip.text = element_text(size = 11, face = "bold"),
        strip.background = element_rect(fill = "gray95", color = NA)
      )
    print(p_genes_clean)
    
    if (length(top_genes_9) > 0) {
      cds_subset <- cds[top_genes_9, ]
      
      p_pseudo_genes <- plot_genes_in_pseudotime(
        cds_subset,
        color_cells_by = "seurat_clusters",
        min_expr = 0.1,
        cell_size = 0.2,
        nrow = 3, ncol = 3) +
        theme(
          legend.position = "bottom",
          strip.text = element_text(size = 10, face = "bold"),
          plot.title = element_text(size = 14, hjust = 0.5, face = "bold")
        ) +
        guides(color = guide_legend(nrow = 2, override.aes = list(size = 3, alpha = 1))) +
        ggtitle("Gene Expression along Pseudotime (by Cluster)")
      print(p_pseudo_genes)
      
      if (has_celltype) {
        p_pseudo_genes_celltype <- plot_genes_in_pseudotime(
          cds_subset,
          color_cells_by = "cell_type",
          min_expr = 0.1,
          cell_size = 0.2,
          nrow = 3, ncol = 3) +
          theme(
            legend.position = "bottom",
            strip.text = element_text(size = 10, face = "bold"),
            plot.title = element_text(size = 14, hjust = 0.5, face = "bold")
          ) +
          guides(color = guide_legend(nrow = 2, override.aes = list(size = 3, alpha = 1))) +
          ggtitle("Gene Expression along Pseudotime (by Cell Type)")
        print(p_pseudo_genes_celltype)
      }
    }
    
    if (nrow(pr_deg) >= 20) {
      top_genes_20 <- head(rownames(pr_deg), 20)
      cds_top20 <- cds[top_genes_20, ]
      
      pt <- SummarizedExperiment::colData(cds)$monocle3_pseudotime
      pt_order <- order(pt)
      
      expr_mat <- as.matrix(SingleCellExperiment::counts(cds_top20)[, pt_order])
      
      expr_scaled <- t(scale(t(log1p(expr_mat))))
      
      n_bins <- 100
      pt_sorted <- pt[pt_order]
      bin_idx <- cut(seq_along(pt_sorted), breaks = n_bins, labels = FALSE)
      
      expr_binned <- sapply(1:n_bins, function(i) {
        rowMeans(expr_scaled[, bin_idx == i, drop = FALSE], na.rm = TRUE)
      })
      
      expr_binned[expr_binned > 2] <- 2
      expr_binned[expr_binned < -2] <- -2
      
      pheatmap(expr_binned,
               cluster_cols = FALSE,
               cluster_rows = TRUE,
               show_colnames = FALSE,
               color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
               main = "Top 20 Trajectory Genes Expression\n(ordered by pseudotime)",
               fontsize_row = 9,
               border_color = NA)
    }
    
    dev.off()
    log_msg(paste("Trajectory dynamics plot saved to:", file.path(output_plots_dir, "gene_dynamics.pdf")))
    
    save_png(p_genes_clean, 
             file.path(output_plots_dir, "top_genes_expression.png"),
             width = 3500, height = 3000)
  }
  
  if (nrow(pr_deg) >= 10) {
    tryCatch({
      log_msg("Finding co-regulated gene modules...")
      
      pr_deg_ids <- head(rownames(pr_deg), min(500, nrow(pr_deg)))
      
      if ("cell_type" %in% colnames(SummarizedExperiment::colData(cds))) {
        cell_groups <- data.frame(
          cell = colnames(cds),
          cell_group = SummarizedExperiment::colData(cds)$cell_type
        )
      } else {
        cell_groups <- data.frame(
          cell = colnames(cds),
          cell_group = SummarizedExperiment::colData(cds)$seurat_clusters
        )
      }
      
      if (length(unique(cell_groups$cell_group)) < 2) {
        log_msg("Cell groups have fewer than 2 classes, skipping gene module analysis", level = "WARN")
      } else {
        gene_modules <- find_gene_modules(cds[pr_deg_ids, ],
                                          resolution = 1e-3)
        
        n_modules <- length(unique(gene_modules$module))
        log_msg(paste("Found", n_modules, "gene modules"))
        
        agg_mat <- aggregate_gene_expression(cds,
                                             gene_group_df = gene_modules,
                                             cell_group_df = cell_groups)
        
        pdf(file.path(output_plots_dir, "gene_modules_heatmap.pdf"), width = 12, height = 10)
        
        row.names(agg_mat) <- paste0("Module_", row.names(agg_mat))
        
        pheatmap(agg_mat,
                 scale = "column",
                 clustering_method = "ward.D2",
                 main = "Gene Module Expression Pattern",
                 fontsize_row = 10,
                 fontsize_col = 9,
                 color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                 border_color = NA,
                 angle_col = 45)
        
        pheatmap(agg_mat,
                 scale = "row",
                 clustering_method = "ward.D2",
                 main = "Gene Module Expression Pattern (Row Scaled)",
                 fontsize_row = 10,
                 fontsize_col = 9,
                 color = colorRampPalette(c("purple4", "white", "darkorange"))(100),
                 border_color = NA,
                 angle_col = 45)
        
        dev.off()
        log_msg(paste("Gene module heatmap saved to:", file.path(output_plots_dir, "gene_modules_heatmap.pdf")))
        
        write.csv(gene_modules,
                  file.path(dirname(output_genes), "gene_modules.csv"),
                  row.names = FALSE)
        log_msg("Gene module information saved")
        
        module_stats <- gene_modules %>%
          group_by(module) %>%
          summarise(
            n_genes = n(),
            genes = paste(head(id, 5), collapse = ", "),
            .groups = "drop"
          ) %>%
          arrange(desc(n_genes))
        
        write.csv(module_stats,
                  file.path(dirname(output_genes), "gene_modules_summary.csv"),
                  row.names = FALSE)
      }
    }, error = function(e) {
      log_msg(paste("Gene module analysis failed, skipping:", e$message), level = "WARN")
    })
  }
  
  log_msg("Saving analysis results...")
  
  saveRDS(cds, output_cds)
  log_msg(paste("CDS object saved to:", basename(output_cds)))
  
  trajectory_summary <- list(
    n_cells = ncol(cds),
    n_genes = nrow(cds),
    n_partitions = length(unique(SummarizedExperiment::colData(cds)$partition)),
    n_clusters = length(unique(SummarizedExperiment::colData(cds)$seurat_clusters)),
    pseudotime_range = range(SummarizedExperiment::colData(cds)$monocle3_pseudotime, na.rm = TRUE),
    n_trajectory_genes = nrow(pr_deg),
    root_method = ifelse(!is.null(root_cells) && length(root_cells) > 0, "Manual", "Auto"),
    analysis_params = list(
      num_dim = num_dim,
      resolution = resolution,
      k = k,
      q_value_threshold = q_value_threshold
    )
  )
  
  saveRDS(trajectory_summary, output_trajectory)
  log_msg(paste("Trajectory analysis summary saved to:", basename(output_trajectory)))
  
  summary_text <- paste(
    "\n========== MONOCLE3 TRAJECTORY ANALYSIS SUMMARY ==========",
    paste("Number of cells:", trajectory_summary$n_cells),
    paste("Number of genes:", trajectory_summary$n_genes),
    paste("Number of partitions:", trajectory_summary$n_partitions),
    paste("Number of clusters:", trajectory_summary$n_clusters),
    paste("Pseudotime range:",
          round(trajectory_summary$pseudotime_range[1], 3), "-",
          round(trajectory_summary$pseudotime_range[2], 3)),
    paste("Number of trajectory genes:", trajectory_summary$n_trajectory_genes),
    paste("Root selection method:", trajectory_summary$root_method),
    paste("Analysis parameters:"),
    paste("  - num_dim:", num_dim),
    paste("  - resolution:", resolution),
    paste("  - k:", k),
    paste("  - q_value_threshold:", q_value_threshold),
    "============================================================\n",
    sep = "\n"
  )
  
  cat(summary_text)
  writeLines(summary_text, file.path(dirname(output_trajectory), "trajectory_summary.txt"))
  
  log_msg("Monocle3 pseudo-time analysis process successfully completed!")
  
}, error = function(e) {
  log_msg(paste("Error:", e$message), level = "ERROR")
  stop(e)
})
