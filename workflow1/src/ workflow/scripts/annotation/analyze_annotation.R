# scripts/annotation_scripts/analyze_annotation.R

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
})
options(future.globals.maxSize = 65 * 1024^3) 
input_file        <- snakemake@input[["auto_annotated"]]
output_seurat     <- snakemake@output[["final_annotated_seurat"]]
output_umap       <- snakemake@output[["umap_annotated"]]
output_composition_plot <- snakemake@output[["composition_plot"]]
output_composition_table <- snakemake@output[["composition_table"]]
output_heatmap    <- snakemake@output[["marker_heatmap"]]
logfc_threshold   <- snakemake@params[["logfc_threshold"]]
min_pct           <- snakemake@params[["min_pct"]]
log_file          <- snakemake@log[[1]]

if (!dir.exists(dirname(log_file))) dir.create(dirname(log_file), recursive = TRUE)
if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] %s", timestamp, level, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  message(line)
}

tryCatch({
  log_msg(paste(rep("=", 70), collapse = ""))
  log_msg("Start the annotation analysis and visualization process (including inter-group comparison)")
  log_msg(paste(rep("=", 70), collapse = ""))
  
  log_msg("Step 1: Load the Seurat object with automatic annotations...")
  seurat_obj <- readRDS(input_file)
  
  log_msg(sprintf("  - Cell count: %d", ncol(seurat_obj)))
  log_msg(sprintf("  - Gene count: %d", nrow(seurat_obj)))
  log_msg(sprintf("  - Cluster count: %d", length(unique(seurat_obj$seurat_clusters))))
  
  required_fields <- c("celltype_main", "celltype_fine")
  missing_fields <- setdiff(required_fields, colnames(seurat_obj@meta.data))
  
  if (length(missing_fields) > 0) {
    log_msg(paste("ERROR: Missing required fields:", paste(missing_fields, collapse = ", ")), level = "ERROR")
    
    if ("singler_labels" %in% colnames(seurat_obj@meta.data)) {
      log_msg("  Attempting to use singler_labels as celltype_fine", level = "WARN")
      if (!"celltype_fine" %in% colnames(seurat_obj@meta.data)) {
        seurat_obj$celltype_fine <- seurat_obj$singler_labels
      }
      if (!"celltype_main" %in% colnames(seurat_obj@meta.data)) {
        seurat_obj$celltype_main <- seurat_obj$singler_labels
        log_msg("  WARNING: celltype_main 未找到，使用 singler_labels 替代", level = "WARN")
      }
    } else {
      stop(paste("Missing required fields and no fallback found:", paste(missing_fields, collapse = ", ")))
    }
  }
  
  if (!"group" %in% colnames(seurat_obj@meta.data)) {
    log_msg("  ERROR: The 'group' column was not detected, unable to perform inter-group comparison", level = "ERROR")
    stop("The 'group' column was not detected")
  }
  
  groups <- if (is.factor(seurat_obj$group)) {
    levels(seurat_obj$group)
  } else {
    sort(unique(seurat_obj$group))
  }
  groups <- groups[!is.na(groups)]
  
  log_msg(sprintf("  - Group count: %d", length(groups)))
  log_msg("  Cell counts per group:")
  group_counts <- table(seurat_obj$group)
  for (g in names(group_counts)) {
    log_msg(sprintf("    - %s: %d cells", g, group_counts[g]))
  }
  
  log_msg(sprintf("  - Main cell type count: %d", length(unique(seurat_obj$celltype_main))))
  log_msg(sprintf("  - Fine cell type count: %d", length(unique(seurat_obj$celltype_fine))))
  
  log_msg("\nStep 2: Find marker genes for each fine cell type...")
  
  DefaultAssay(seurat_obj) <- "SCT"
  
  log_msg("  Running PrepSCTFindMarkers()...")
  seurat_obj <- tryCatch({
    PrepSCTFindMarkers(seurat_obj, verbose = FALSE)
  }, error = function(e) {
    log_msg(paste("  PrepSCTFindMarkers failed, continuing with existing SCT data:", e$message), level = "WARN")
    seurat_obj
  })
  
  log_msg("  Setting cell identity to celltype_fine...")
  Idents(seurat_obj) <- "celltype_fine"
  
  n_celltypes <- length(unique(seurat_obj$celltype_fine))
  if (n_celltypes < 2) {
    log_msg("  WARNING: Only one cell type found, skipping marker gene search", level = "WARN")
    all_markers <- data.frame()
    top10_markers <- data.frame()
  } else {
    log_msg(sprintf("  Running FindAllMarkers() (total %d cell types)...", n_celltypes))

    all_markers <- tryCatch({
      FindAllMarkers(
        seurat_obj,
        assay = "SCT",
        only.pos = TRUE,
        min.pct = min_pct,
        logfc.threshold = logfc_threshold,
        verbose = FALSE
      )
    }, error = function(e) {
      log_msg(paste("  ERROR in FindAllMarkers:", e$message), level = "ERROR")
      return(data.frame())
    })
    
    if (nrow(all_markers) > 0) {
      log_msg(sprintf("  Found %d marker genes", nrow(all_markers)))
      
      top10_markers <- all_markers %>%
        group_by(cluster) %>%
        dplyr::slice_max(n = 10, order_by = avg_log2FC)
      
      log_msg(sprintf("  Extracted Top 10 markers for each cell type"))
      log_msg("  Marker gene data will be saved to the Seurat object (CSV file not generated)")
    } else {
      log_msg("  WARNING: No marker genes found", level = "WARN")
      top10_markers <- data.frame()
    }
  }
  
  log_msg("\nStep 3: Prepare visualization color schemes...")
  
  saved_idents <- Idents(seurat_obj)
  
  Idents(seurat_obj) <- "seurat_clusters"
  cluster_levels <- levels(Idents(seurat_obj))
  if (is.null(cluster_levels)) {
    cluster_levels <- sort(unique(as.character(seurat_obj$seurat_clusters)))
  }
  cluster_cols <- setNames(hue_pal()(length(cluster_levels)), cluster_levels)
  
  main_levels <- sort(unique(seurat_obj$celltype_main))
  main_cols <- setNames(hue_pal(h = c(15, 375), l = 65, c = 100)(length(main_levels)), 
                        main_levels)
  
  fine_levels <- sort(unique(seurat_obj$celltype_fine))
  fine_cols <- setNames(hue_pal(h = c(0, 360), l = 70, c = 90)(length(fine_levels)), 
                        fine_levels)
  
  group_cols <- setNames(brewer.pal(max(3, length(groups)), "Set2")[1:length(groups)], 
                         groups)
  
  log_msg(sprintf("  - Color schemes initialized: %d clusters, %d main types, %d fine types, %d groups", 
                  length(cluster_levels), 
                  length(main_levels), 
                  length(fine_levels),
                  length(groups)))
  
  log_msg("\nStep 4: Generate marker gene heatmap...")
  
  sapply(c(output_seurat, output_umap, output_composition_plot, 
           output_composition_table, output_heatmap), 
         function(f) {
           d <- dirname(f)
           if (!dir.exists(d)) dir.create(d, recursive = TRUE)
         })
  
  if (nrow(top10_markers) > 0) {
    max_cells_per_type <- 100
    total_celltypes <- length(unique(seurat_obj$celltype_fine))
    
    if (ncol(seurat_obj) > max_cells_per_type * total_celltypes) {
      log_msg(sprintf("  Object is large, showing at most %d cells per cell type", max_cells_per_type))
      
      cells_to_plot <- unlist(lapply(
        split(rownames(seurat_obj@meta.data), seurat_obj@meta.data$celltype_fine), 
        function(cells) {
          n_sample <- min(max_cells_per_type, length(cells))
          sample(cells, n_sample)
        }
      ))
      
      seurat_subset <- subset(seurat_obj, cells = cells_to_plot)
      log_msg(sprintf("  Sampled cell count: %d", ncol(seurat_subset)))
    } else {
      seurat_subset <- seurat_obj
    }
    
    Idents(seurat_subset) <- "celltype_fine"
    
    log_msg("  Plotting heatmap...")
    pdf(output_heatmap, width = 16, height = 20)
    
    tryCatch({
      p_heatmap <- DoHeatmap(
        seurat_subset, 
        features = top10_markers$gene,
        size = 3,
        angle = 0,
        hjust = 0.5
      ) + 
        NoLegend() +
        theme(axis.text.y = element_text(size = 6)) +
        labs(title = "Top 10 Marker Genes per Cell Type")
      
      print(p_heatmap)
    }, error = function(e) {
      log_msg(paste("  WARNING: Heatmap generation failed:", e$message), level = "WARN")
      plot.new()
      text(0.5, 0.5, paste("Heatmap generation failed:\n", e$message), cex = 1.5)
    })
    
    dev.off()
    log_msg(sprintf("  Heatmap saved: %s", basename(output_heatmap)))
  } else {
    log_msg("  Skipping heatmap generation (no marker genes)", level = "WARN")
    pdf(output_heatmap, width = 10, height = 6)
    plot.new()
    text(0.5, 0.5, "No available marker genes, heatmap not generated", cex = 1.2)
    dev.off()
  }
  log_msg("\nStep 5: Calculate cell type composition data...")
  
  sample_col <- if ("sample_id" %in% colnames(seurat_obj@meta.data)) {
    "sample_id"
  } else if ("orig.ident" %in% colnames(seurat_obj@meta.data)) {
    "orig.ident"
  } else {
    log_msg("  WARNING: Sample ID column not found, using group as sample", level = "WARN")
    "group"
  }
  
  log_msg(sprintf("  Using sample column: %s", sample_col))
  
  composition_main_sample <- as.data.frame(
    prop.table(
      table(seurat_obj$celltype_main, seurat_obj@meta.data[[sample_col]]), 
      margin = 2
    )
  )
  colnames(composition_main_sample) <- c("CellType", "Sample", "Proportion")
  composition_main_sample$Level <- "Major"
  
  composition_fine_sample <- as.data.frame(
    prop.table(
      table(seurat_obj$celltype_fine, seurat_obj@meta.data[[sample_col]]), 
      margin = 2
    )
  )
  colnames(composition_fine_sample) <- c("CellType", "Sample", "Proportion")
  composition_fine_sample$Level <- "Fine"
  
  composition_main_group <- as.data.frame(
    prop.table(
      table(seurat_obj$celltype_main, seurat_obj$group), 
      margin = 2
    )
  )
  colnames(composition_main_group) <- c("CellType", "Group", "Proportion")
  
  composition_fine_group <- as.data.frame(
    prop.table(
      table(seurat_obj$celltype_fine, seurat_obj$group), 
      margin = 2
    )
  )
  colnames(composition_fine_group) <- c("CellType", "Group", "Proportion")
  
  cellcount_main_group <- as.data.frame(
    table(seurat_obj$celltype_main, seurat_obj$group)
  )
  colnames(cellcount_main_group) <- c("CellType", "Group", "Count")
  
  cellcount_fine_group <- as.data.frame(
    table(seurat_obj$celltype_fine, seurat_obj$group)
  )
  colnames(cellcount_fine_group) <- c("CellType", "Group", "Count")
  
  composition_data <- rbind(composition_main_sample, composition_fine_sample)
  
  write.csv(composition_data, file = output_composition_table, row.names = FALSE)
  log_msg(sprintf("  Composition table saved: %s", basename(output_composition_table)))
  
  group_composition_file <- sub("\\.csv$", "_by_group.csv", output_composition_table)
  composition_group_all <- rbind(
    cbind(composition_main_group, Level = "Major"),
    cbind(composition_fine_group, Level = "Fine")
  )
  write.csv(composition_group_all, file = group_composition_file, row.names = FALSE)
  log_msg(sprintf("  Group composition table saved: %s", basename(group_composition_file)))
  
  log_msg("\nStep 6: Generate integrated visualization plots (UMAP + inter-group comparisons, all in one file)...")
  
  pdf(output_umap, width = 18, height = 10, onefile = TRUE)
  
  Idents(seurat_obj) <- "seurat_clusters"
  p_umap_clusters <- DimPlot(
    seurat_obj, 
    reduction = "umap", 
    label = TRUE, 
    repel = TRUE, 
    cols = cluster_cols,
    label.size = 4,
    pt.size = 0.5
  ) +
    labs(title = "Original Clusters") + 
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  Idents(seurat_obj) <- "celltype_main"
  p_umap_main <- DimPlot(
    seurat_obj, 
    reduction = "umap", 
    label = TRUE, 
    repel = TRUE, 
    cols = main_cols,
    label.size = 5,
    pt.size = 0.5
  ) +
    labs(title = "Major Cell Types") + 
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  Idents(seurat_obj) <- "celltype_fine"
  p_umap_fine <- DimPlot(
    seurat_obj, 
    reduction = "umap", 
    label = TRUE, 
    repel = TRUE, 
    cols = fine_cols,
    label.size = 4,
    pt.size = 0.5
  ) +
    labs(title = "Cell Subtypes") + 
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  print(p_umap_clusters | p_umap_main | p_umap_fine)
  
  p_umap_group <- DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "group",
    cols = group_cols,
    pt.size = 0.5
  ) +
    labs(title = "UMAP Colored by Group") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          legend.position = "right")
  
  print(p_umap_group)
  
  Idents(seurat_obj) <- "celltype_main"
  p_umap_split_main <- DimPlot(
    seurat_obj,
    reduction = "umap",
    split.by = "group",
    label = TRUE,
    repel = TRUE,
    cols = main_cols,
    label.size = 4,
    pt.size = 0.3,
    ncol = length(groups)
  ) +
    labs(title = "Major Cell Types Split by Group") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  print(p_umap_split_main)
  
  Idents(seurat_obj) <- "celltype_fine"
  p_umap_split_fine <- DimPlot(
    seurat_obj,
    reduction = "umap",
    split.by = "group",
    label = TRUE,
    repel = TRUE,
    cols = fine_cols,
    label.size = 3,
    pt.size = 0.3,
    ncol = length(groups)
  ) +
    labs(title = "Cell Subtypes Split by Group") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  print(p_umap_split_fine)
  
  log_msg(sprintf("  Page 5+: Generating UMAP for %d groups...", length(groups)))
  
  for (g in groups) {
    log_msg(sprintf("    - Plotting group: %s", g))
    sub_obj <- subset(seurat_obj, subset = group == g)
    
    if (ncol(sub_obj) == 0) {
      log_msg(sprintf("      Group '%s' has no cells, skipping", g), level = "WARN")
      next
    }
    
    Idents(sub_obj) <- "celltype_main"
    p_g_main <- DimPlot(
      sub_obj, 
      reduction = "umap", 
      label = TRUE,
      repel = TRUE, 
      cols = main_cols, 
      label.size = 5,
      pt.size = 0.8
    ) +
      NoLegend() + 
      labs(title = paste0("Major Types - Group: ", g)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    
    Idents(sub_obj) <- "celltype_fine"
    p_g_fine <- DimPlot(
      sub_obj, 
      reduction = "umap", 
      label = TRUE,
      repel = TRUE, 
      cols = fine_cols, 
      label.size = 4,
      pt.size = 0.8
    ) +
      NoLegend() + 
      labs(title = paste0("Subtypes - Group: ", g)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    
    print(p_g_main | p_g_fine)
  }
  
  log_msg("  Generating inter-group comparison plots...")
  
  log_msg("  Inter-group comparison plot 1: Absolute cell counts of major cell types...")
  p_count_main <- ggplot(cellcount_main_group, 
                         aes(x = CellType, y = Count, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_manual(values = group_cols) +
    labs(
      title = "Absolute Cell Counts of Major Cell Types by Group", 
      x = "Cell Type", 
      y = "Cell Count",
      fill = "Group"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    )
  
  print(p_count_main)
  
  log_msg("  Inter-group comparison plot 2: Top 20 cell subtypes absolute cell counts...")
  top20_subtypes <- cellcount_fine_group %>%
    group_by(CellType) %>%
    summarise(total = sum(Count), .groups = 'drop') %>%
    arrange(desc(total)) %>%
    head(20) %>%
    pull(CellType)
  
  cellcount_fine_top20 <- cellcount_fine_group %>%
    filter(CellType %in% top20_subtypes)
  
  p_count_fine <- ggplot(cellcount_fine_top20, 
                         aes(x = reorder(CellType, -Count), y = Count, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_manual(values = group_cols) +
    labs(
      title = "Absolute Cell Counts of Top 20 Cell Subtypes by Group", 
      x = "Cell Subtype", 
      y = "Cell Count",
      fill = "Group"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "top"
    )
  
  print(p_count_fine)
  
  log_msg("  Inter-group comparison plot 3: Major cell type inter-group proportion comparison...")
  p_comp_main_group <- ggplot(composition_main_group, 
                              aes(x = Group, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill", width = 0.6) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    scale_fill_manual(values = main_cols) +
    labs(
      title = "Major Cell Type Composition Comparison Between Groups", 
      x = "Group", 
      y = "Proportion",
      fill = "Cell Type"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
  
  print(p_comp_main_group)
  
  log_msg("  Inter-group comparison plot 4: Cell subtype inter-group proportion comparison...")
  p_comp_fine_group <- ggplot(composition_fine_group, 
                              aes(x = Group, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill", width = 0.6) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    scale_fill_manual(values = fine_cols) +
    labs(
      title = "Cell Subtype Composition Comparison Between Groups", 
      x = "Group", 
      y = "Proportion",
      fill = "Cell Subtype"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.text = element_text(size = 7)
    )
  
  print(p_comp_fine_group)
  
  log_msg("  Inter-group comparison plot 5: Major cell type proportions by group...")
  p_comp_main_facet <- ggplot(composition_main_group, 
                              aes(x = CellType, y = Proportion, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = group_cols) +
    labs(
      title = "Major Cell Type Proportions by Group", 
      x = "Cell Type", 
      y = "Proportion",
      fill = "Group"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    )
  
  print(p_comp_main_facet)
  
  log_msg("  Inter-group comparison plot 6: Top 15 cell subtypes proportions by group...")
  top15_subtypes <- composition_fine_group %>%
    group_by(CellType) %>%
    summarise(total = sum(Proportion), .groups = 'drop') %>%
    arrange(desc(total)) %>%
    head(15) %>%
    pull(CellType)
  
  composition_fine_top15 <- composition_fine_group %>%
    filter(CellType %in% top15_subtypes)
  
  p_comp_fine_facet <- ggplot(composition_fine_top15, 
                              aes(x = reorder(CellType, -Proportion), y = Proportion, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = group_cols) +
    labs(
      title = "Top 15 Cell Subtype Proportions by Group", 
      x = "Cell Subtype", 
      y = "Proportion",
      fill = "Group"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "top"
    )
  
  print(p_comp_fine_facet)
  
  log_msg("  Inter-group comparison plot 7: Major cell type composition heatmap...")
  
  heatmap_main_data <- composition_main_group %>%
    pivot_wider(names_from = Group, values_from = Proportion, values_fill = 0) %>%
    as.data.frame()
  
  rownames(heatmap_main_data) <- heatmap_main_data$CellType
  heatmap_main_data$CellType <- NULL
  heatmap_main_mat <- as.matrix(heatmap_main_data)
  if (length(unique(as.vector(heatmap_main_mat))) < 2) {
    heatmap_main_mat[1, 1] <- heatmap_main_mat[1, 1] + 1e-06
  }
  
  pheatmap(
    heatmap_main_mat * 100,
    cluster_rows = nrow(heatmap_main_mat) >= 2,
    cluster_cols = ncol(heatmap_main_mat) >= 2,
    color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
    display_numbers = TRUE,
    number_format = "%.1f",
    main = "Major Cell Type Composition Heatmap (%)",
    fontsize = 12,
    fontsize_number = 10,
    angle_col = 0
  )
  
  log_msg("  Inter-group comparison plot 8: Heatmap of the 20 most variable cell subtypes...")
  
  subtype_variance <- composition_fine_group %>%
    group_by(CellType) %>%
    summarise(
      mean_prop = mean(Proportion),
      sd_prop = sd(Proportion),
      cv = ifelse(mean_prop > 0, sd_prop / mean_prop, 0),
      .groups = 'drop'
    ) %>%
    arrange(desc(cv)) %>%
    head(20)
  
  heatmap_fine_data <- composition_fine_group %>%
    filter(CellType %in% subtype_variance$CellType) %>%
    pivot_wider(names_from = Group, values_from = Proportion, values_fill = 0) %>%
    as.data.frame()
  
  rownames(heatmap_fine_data) <- heatmap_fine_data$CellType
  heatmap_fine_data$CellType <- NULL
  heatmap_fine_mat <- as.matrix(heatmap_fine_data)
  if (length(unique(as.vector(heatmap_fine_mat))) < 2) {
    heatmap_fine_mat[1, 1] <- heatmap_fine_mat[1, 1] + 1e-06
  }
  
  pheatmap(
    heatmap_fine_mat * 100,
    cluster_rows = nrow(heatmap_fine_mat) >= 2,
    cluster_cols = ncol(heatmap_fine_mat) >= 2,
    color = colorRampPalette(c("white", "orange", "red"))(100),
    display_numbers = TRUE,
    number_format = "%.1f",
    main = "Top 20 Variable Cell Subtype Composition Heatmap (%)",
    fontsize = 11,
    fontsize_number = 9,
    angle_col = 0
  )
  
  if (length(groups) == 2) {
    log_msg("  Inter-group comparison plot 9: Cell type differences between two groups...")
    
    composition_main_wide <- composition_main_group %>%
      pivot_wider(names_from = Group, values_from = Proportion, values_fill = 0)
    
    composition_main_wide$Difference <- composition_main_wide[[groups[1]]] - 
      composition_main_wide[[groups[2]]]
    
    composition_main_wide <- composition_main_wide %>%
      arrange(desc(abs(Difference)))
    
    p_diff <- ggplot(composition_main_wide, 
                     aes(x = reorder(CellType, Difference), y = Difference * 100)) +
      geom_bar(stat = "identity", aes(fill = Difference > 0), width = 0.7) +
      scale_fill_manual(values = c("steelblue", "coral"), 
                        labels = c(paste("Higher in", groups[2]), 
                                   paste("Higher in", groups[1]))) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
      labs(
        title = sprintf("Cell Type Proportion Differences (%s vs %s)", 
                        groups[1], groups[2]),
        x = "Cell Type",
        y = "Proportion Difference (%)",
        fill = "Enriched in"
      ) +
      coord_flip() +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.y = element_text(size = 11),
        legend.position = "top",
        legend.title = element_text(face = "bold")
      )
    
    print(p_diff)
  }
  
  if (sample_col != "group" && length(unique(seurat_obj@meta.data[[sample_col]])) > length(groups)) {
    log_msg("  Inter-group comparison plot 10: Cell type distribution at the sample level...")
    
    sample_celltype_prop <- seurat_obj@meta.data %>%
      group_by(across(all_of(c(sample_col, "group", "celltype_main")))) %>%
      summarise(n = n(), .groups = 'drop') %>%
      group_by(across(all_of(sample_col))) %>%
      mutate(proportion = n / sum(n))
    
    major_celltypes <- names(sort(table(seurat_obj$celltype_main), decreasing = TRUE)[1:min(6, length(main_levels))])
    
    sample_celltype_prop_filtered <- sample_celltype_prop %>%
      filter(celltype_main %in% major_celltypes)
    
    p_violin <- ggplot(sample_celltype_prop_filtered, 
                       aes(x = group, y = proportion * 100, fill = group)) +
      geom_violin(alpha = 0.6) +
      geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
      geom_jitter(width = 0.1, size = 1.5, alpha = 0.5) +
      facet_wrap(~ celltype_main, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = group_cols) +
      labs(
        title = "Distribution of Major Cell Type Proportions Across Samples by Group",
        x = "Group",
        y = "Proportion (%)",
        fill = "Group"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "lightgray"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none"
      )
    
    print(p_violin)
  }
  
  dev.off()  
  log_msg("\nGenerate a composition visualization chart...")
  
  pdf(output_composition_plot, width = 14, height = 10)
  
  p_comp_main_sample <- ggplot(composition_main_sample, 
                               aes(x = Sample, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    scale_fill_manual(values = main_cols) +
    labs(
      title = "Major Cell Type Composition by Sample", 
      x = "Sample", 
      y = "Proportion",
      fill = "Cell Type"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    )
  
  p_comp_fine_sample <- ggplot(composition_fine_sample, 
                               aes(x = Sample, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    scale_fill_manual(values = fine_cols) +
    labs(
      title = "Cell Subtype Composition by Sample", 
      x = "Sample", 
      y = "Proportion",
      fill = "Cell Subtype"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.text = element_text(size = 8)
    )
  
  print(p_comp_main_sample)
  print(p_comp_fine_sample)
  
  dev.off()
  log_msg(sprintf("  Composition plot saved: %s", basename(output_composition_plot)))
  
  log_msg("\nStep 8: Save the final Seurat object...")
  
  seurat_obj@misc$top10_markers <- top10_markers
  seurat_obj@misc$all_markers <- all_markers
  seurat_obj@misc$composition_data <- composition_data
  seurat_obj@misc$composition_by_group <- composition_group_all
  seurat_obj@misc$cellcount_by_group <- list(
    main = cellcount_main_group,
    fine = cellcount_fine_group
  )
  seurat_obj@misc$analysis_params <- list(
    logfc_threshold = logfc_threshold,
    min_pct = min_pct,
    analysis_date = Sys.time()
  )
  
  Idents(seurat_obj) <- saved_idents
  
  saveRDS(seurat_obj, file = output_seurat, compress = FALSE)
  log_msg(sprintf("  Final object saved: %s", basename(output_seurat)))
  log_msg(sprintf("  Object size: %s", format(object.size(seurat_obj), units = "MB")))
  
  log_msg("\n" %>% paste(rep("=", 70), collapse = ""))
  log_msg("Annotation analysis and visualization workflow completed successfully!")
  log_msg(paste(rep("=", 70), collapse = ""))
  
  log_msg("\nSummary:")
  log_msg(sprintf("  - Total cells: %d", ncol(seurat_obj)))
  log_msg(sprintf("  - Total clusters: %d", length(unique(seurat_obj$seurat_clusters))))
  log_msg(sprintf("  - Number of groups: %d", length(groups)))
  log_msg(sprintf("  - Number of major cell types: %d", length(unique(seurat_obj$celltype_main))))
  log_msg(sprintf("  - Number of cell subtypes: %d", length(unique(seurat_obj$celltype_fine))))
  log_msg(sprintf("  - Total marker genes: %d", nrow(all_markers)))
  
  log_msg("\nOutput files:")
  log_msg(sprintf("  1. Final Seurat object: %s", output_seurat))
  log_msg(sprintf("  2. Integrated visualization (UMAP + inter-group comparisons): %s", output_umap))
  log_msg(sprintf("  3. Marker heatmap: %s", output_heatmap))
  log_msg(sprintf("  4. Composition plot: %s", output_composition_plot))
  log_msg(sprintf("  5. Composition table (sample-level): %s", output_composition_table))
  log_msg(sprintf("  6. Composition table (group-level): %s", group_composition_file))
  
  log_msg("\nNotes:")
  log_msg("  - Marker gene data has been saved to the Seurat object's misc slot (not output as CSV)")
  log_msg("  - All marker genes can be accessed via seurat_obj@misc$all_markers")
  log_msg("  - Top 10 markers can be accessed via seurat_obj@misc$top10_markers")
  log_msg("  - Composition data is saved to both the CSV file and the Seurat object")
  
  log_msg("\nIntegrated visualization PDF content:")
  log_msg("  - UMAP section: Cluster/main type/subtype comparison, group coloring, Split view, detailed views for each group")
  log_msg("  - Inter-group comparisons: Absolute numbers, proportions, heatmaps, differential analysis, distribution violin plots")
  
  log_msg("\n" %>% paste(rep("=", 70), collapse = ""))
  
}, error = function(e) {
  log_msg(paste(rep("=", 70), collapse = ""), level = "ERROR")
  log_msg(paste("Fatal error:", e$message), level = "ERROR")
  log_msg(paste(rep("=", 70), collapse = ""), level = "ERROR")
  
  log_msg("\nError traceback:", level = "ERROR")
  traceback_info <- capture.output(traceback())
  for (line in traceback_info) {
    log_msg(line, level = "ERROR")
  }
  
  stop(e)
})
