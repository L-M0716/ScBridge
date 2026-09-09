# run_auto_annotation.R
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(UCell)
  library(BiocParallel)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(data.table)
  library(stringr)
  library(tibble)
})

input_file     <- snakemake@input[["clustered_seurat"]]
panglaodb_file <- snakemake@input[["markerlist"]]
output_seurat  <- snakemake@output[["auto_annotated"]]
output_report  <- snakemake@output[["annotation_report"]]
output_scores  <- snakemake@output[["confidence_scores"]]
log_file       <- snakemake@log[[1]]
organism       <- snakemake@params[["organism"]]

# ========== parameters ==========
UCELL_MAXRANK    <- 1500
UCELL_MINMARKERS <- 2
UCELL_MARGIN     <- 0.03   # Intra-subtype competition: accept only if highest - second highest >= this value

SINGLER_GENES       <- "de"
SINGLER_SD_THRESH   <- 1
SINGLER_DE_METHOD   <- "classic"
SINGLER_QUANTILE    <- 0.8
SINGLER_TUNE_THRESH <- 0.05
SINGLER_ASSAY       <- "logcounts"
PRUNE_N_MADS        <- 3
PRUNE_MIN_DIFF_MED  <- -Inf
PRUNE_MIN_DIFF_NEXT <- 0
HPCA_CONF_THRESH <- 0.20
log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H%M%S")
  line <- sprintf("[%s] [%s] %s", timestamp, level, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  cat(line, "\n")
}
get_gene_converter <- function(organism) {
  if (tolower(organism) %in% c("mouse", "mm", "mus musculus")) {
    return(function(g) sapply(g, function(x) {
      if (is.na(x) || x == "") return(NA_character_)
      paste0(toupper(substr(x,1,1)), tolower(substr(x,2,nchar(x))))
    }, USE.NAMES = FALSE))
  } else {
    return(function(g) sapply(g, function(x) {
      if (is.na(x) || x == "") return(NA_character_)
      toupper(x)
    }, USE.NAMES = FALSE))
  }
}

filter_by_species <- function(db, organism) {
  sp <- tolower(trimws(as.character(db$species)))
  if (tolower(organism) %in% c("human", "hs", "homo sapiens")) {
    keep <- sp %in% c("hs","human","homo sapiens") |
      grepl("hs", sp, fixed=TRUE) | grepl("human", sp, fixed=TRUE) |
      grepl("homo", sp, fixed=TRUE)
  } else {
    keep <- sp %in% c("mm","mouse","mus musculus") |
      grepl("mm", sp, fixed=TRUE) | grepl("mouse", sp, fixed=TRUE) |
      grepl("mus", sp, fixed=TRUE)
  }
  db[keep, ]
}

safe_hue_pal <- function(n, ...) {
  if (n <= 0) return(character(0))
  if (n == 1) return("#F8766D")
  hue_pal(...)(n)
}



# ---- human HPCA Mapping table ----
HPCA_HARDCODE_MAP <- c(
  "dc"                  = "Dendritic cells",
  "gmp"                 = "Progenitors",
  "mep"                 = "Progenitors",
  "cmp"                 = "Progenitors",
  "hsc  g csf"          = "Progenitors",
  "hsc cd34+"           = "Progenitors",
  "hsc"                 = "Progenitors",
  "erythroblast"        = "Progenitors",
  "megakaryocytes"      = "Progenitors",
  "platelets"           = "Platelets",
  "myelocyte"           = "Granulocytes",
  "pro myelocyte"       = "Granulocytes",
  "neutrophils"         = "Granulocytes",
  "eosinophils"         = "Granulocytes",
  "mast cells"          = "Basophils",
  "pre b cell cd34 "    = "Progenitors",
  "pro b cell cd34+"    = "Progenitors",
  "bm"                  = "Progenitors",
  "bm  prog."           = "Progenitors",
  "msc"                 = "Stromal cells",
  "chondrocytes"        = "Stromal cells",
  "osteoblasts"         = "Stromal cells",
  "smooth muscle cells" = "Smooth muscle cells",
  "ips cells"           = NA_character_,
  "embryonic stem cells"= NA_character_
)

# ---- mouse MouseRNAseqData Mapping table ----
MOUSE_HARDCODE_MAP <- c(
  "natural killer cells" = "NK cells",
  "nk cell"              = "NK cells",
  "fibroblasts"          = "Stromal cells",
  "msc"                  = "Stromal cells",
  "chondrocytes"         = "Stromal cells",
  "osteoblasts"          = "Stromal cells",
  "erythrocytes"         = "Erythrocytes",
  "smooth muscle cells"  = "Smooth muscle cells",
  "granulocytes"         = "Granulocytes",
  "mast cells"           = "Basophils",
  "ips cells"            = NA_character_,
  "embryonic stem cells" = NA_character_
)

normalize_label <- function(x) {
  x <- tolower(x)
  x <- gsub("[_\\-\\.]+", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  x
}

fuzzy_match_hpca <- function(ref_label, marker_main_types,
                              hardcode_map = HPCA_HARDCODE_MAP) {
  norm_ref     <- normalize_label(ref_label)
  norm_markers <- sapply(marker_main_types, normalize_label)

  if (norm_ref %in% names(hardcode_map)) {
    target <- hardcode_map[[norm_ref]]
    if (is.na(target)) return(NA_character_)
    if (target %in% marker_main_types) return(target)
  }

  exact <- names(norm_markers)[norm_markers == norm_ref]
  if (length(exact) > 0) return(exact[1])

  norm_ref_ns     <- sub("s$", "", norm_ref)
  norm_markers_ns <- sub("s$", "", norm_markers)
  sing <- names(norm_markers_ns)[norm_markers_ns == norm_ref_ns]
  if (length(sing) > 0) return(sing[1])

  for (orig in names(norm_markers)) {
    m    <- norm_markers[[orig]]
    m_ns <- sub("s$", "", m)
    if (grepl(paste0("^", m_ns), norm_ref) ||
        grepl(paste0("^", norm_ref_ns), m)) {
      return(orig)
    }
  }

  hit <- agrep(norm_ref, norm_markers, max.distance = 0.25, value = FALSE)
  if (length(hit) > 0) return(names(norm_markers)[hit[1]])

  return(NA_character_)
}

tryCatch({

  convert_gene_format <- get_gene_converter(organism)
  is_mouse <- tolower(organism) %in% c("mouse", "mm", "mus musculus", "mus_musculus", "mmu")
  species_label <- if (is_mouse) "Mouse" else "Human"

  ACTIVE_HARDCODE_MAP <- if (is_mouse) MOUSE_HARDCODE_MAP else HPCA_HARDCODE_MAP

  log_msg("=== Hierarchical cell type annotation ===")
  log_msg("[Step 1] Loading Seurat object...")
  seurat_obj <- readRDS(input_file)

  available_assays <- names(seurat_obj@assays)
  log_msg(sprintf("   Assay detected: %s", paste(available_assays, collapse=", ")))

  if ("RNA" %in% available_assays) {
    use_assay <- "RNA"
    log_msg("   assay Select: RNA")
    DefaultAssay(seurat_obj) <- "RNA"

    rna_data_raw <- tryCatch(
      GetAssayData(seurat_obj, assay="RNA", layer="data"),
      error = function(e) NULL
    )
    needs_norm <- TRUE
    if (!is.null(rna_data_raw) && nrow(rna_data_raw) > 0 && ncol(rna_data_raw) > 0) {
      sample_genes <- sample(rownames(rna_data_raw), min(100, nrow(rna_data_raw)))
      sample_cells <- seq_len(min(20, ncol(rna_data_raw)))
      rna_max <- tryCatch(
        max(rna_data_raw[sample_genes, sample_cells], na.rm=TRUE),
        error = function(e) NA_real_
      )
      needs_norm <- !is.finite(rna_max) || rna_max > 20
    }
    if (needs_norm) {
      log_msg("   The RNA data layer is empty or not normalized, executing NormalizeData...", level="WARN")
      seurat_obj <- NormalizeData(seurat_obj, assay="RNA",
                                  normalization.method="LogNormalize",
                                  scale.factor=10000, verbose=FALSE)
    }
    tryCatch({
      seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
      log_msg("   RNA assay JoinLayers completed.")
    }, error = function(e) {
      log_msg(sprintf("   RNA JoinLayers skipped: %s", e$message), "WARN")
    })

  } else if ("SCT" %in% available_assays) {
    use_assay <- "SCT"
    log_msg("   assay Select: SCT", level="WARN")
    DefaultAssay(seurat_obj) <- "SCT"
  } else {
    stop(sprintf("No available assay (detected: %s), RNA or SCT required.",
                 paste(available_assays, collapse=", ")))
  }

  log_msg(sprintf("   Cells %d, Genes %d (assay: %s)",
                  ncol(seurat_obj), nrow(seurat_obj), use_assay))

  original_clusters <- as.character(seurat_obj$seurat_clusters)
  cluster_levels <- as.character(levels(seurat_obj$seurat_clusters))
  if (is.null(cluster_levels) || length(cluster_levels) == 0)
    cluster_levels <- as.character(sort(as.numeric(unique(original_clusters))))
  log_msg(sprintf("   Cluster  %d (%s)", length(cluster_levels),
                  paste(cluster_levels, collapse=", ")))

  log_msg("\n[Step 2] Loading marker database...")
  db <- fread(panglaodb_file, sep="\t", quote="", header=TRUE)
  if (length(colnames(db)) > 0)
    colnames(db)[1] <- sub("^\uFEFF", "", colnames(db)[1])

  if (!("marker_gene" %in% colnames(db)) && "official gene symbol" %in% colnames(db))
    db[, marker_gene := get("official gene symbol")]
  if (!("marker_gene" %in% colnames(db)) && "official_gene_symbol" %in% colnames(db))
    db[, marker_gene := get("official_gene_symbol")]
  if (!("main_cell_type" %in% colnames(db)) && "cell type" %in% colnames(db))
    db[, main_cell_type := get("cell type")]
  if (!("main_cell_type" %in% colnames(db)) && "cell_type" %in% colnames(db))
    db[, main_cell_type := get("cell_type")]
  if (!("sub_cell_type" %in% colnames(db)))
    db[, sub_cell_type := NA_character_]

  db_species <- filter_by_species(db, organism)
  log_msg(sprintf("   Filtered by species '%s': %d rows", organism, nrow(db_species)))
  if (nrow(db_species) == 0) stop(sprintf("No matching records for species '%s' in the database", organism))

  db_species$gene_converted <- convert_gene_format(db_species$marker_gene)
  db_species <- db_species[!is.na(gene_converted) & gene_converted != "", ]
  seurat_genes <- rownames(seurat_obj)
  n_before <- length(unique(db_species$gene_converted))
  db_species <- db_species[gene_converted %in% seurat_genes, ]
  n_after  <- length(unique(db_species$gene_converted))
  log_msg(sprintf("   Database genes %d, matched to object %d (%.1f%%)",
                  n_before, n_after, 100 * n_after / max(1, n_before)))

  main_only_df <- db_species[is.na(sub_cell_type) | trimws(sub_cell_type) == "", ]
  sub_df       <- db_species[!is.na(sub_cell_type) & trimws(sub_cell_type) != "", ]

  marker_main_types <- unique(db_species$main_cell_type)
  log_msg(sprintf("   marker types: %s",
                  paste(marker_main_types, collapse=", ")))

  fine_markers <- list()
  fine_to_main_db <- character(0)
  if (nrow(sub_df) > 0) {
    fine_markers <- split(sub_df$gene_converted, sub_df$sub_cell_type)
    fine_markers <- lapply(fine_markers, unique)
    fine_markers <- fine_markers[sapply(fine_markers, length) >= UCELL_MINMARKERS]
    ft_map <- sub_df %>% distinct(sub_cell_type, main_cell_type)
    ft_map <- ft_map[!duplicated(ft_map$sub_cell_type), ]
    fine_to_main_db <- setNames(ft_map$main_cell_type, ft_map$sub_cell_type)
  }
  log_msg(sprintf("   Subtype definitions: %d (marker >= %d)",
                  length(fine_markers), UCELL_MINMARKERS))

  log_msg("\n========== [Layer 1] SingleR + HPCA → Broad Cell Type Annotation ==========")
  log_msg("Reference dataset: HumanPrimaryCellAtlasData (Mabbott et al. 2013 BMC Genomics)")
  log_msg("Coverage: Immune, Epithelial, Endothelial, Mesenchymal, Neural, Hepatocyte and other major lineages")

  log_msg("[L1-1] changing SingleCellExperiment...")
  sce_test <- as.SingleCellExperiment(seurat_obj, assay = use_assay)
  if (!("logcounts" %in% SummarizedExperiment::assayNames(sce_test))) {
    SummarizedExperiment::assay(sce_test, "logcounts") <-
      SummarizedExperiment::assay(sce_test, "data")
  }

  log_msg("[L1-2] Loading reference dataset...")
  if (!is_mouse) {
    ref_broad <- celldex::HumanPrimaryCellAtlasData(ensembl = FALSE)
    log_msg("   Human species: HumanPrimaryCellAtlasData (Mabbott et al. 2013 BMC Genomics)")
    log_msg("   Coverage: Immune, Epithelial, Endothelial, Mesenchymal, Neural, Hepatocyte and other major lineages")
  } else {
    ref_broad <- celldex::MouseRNAseqData(ensembl = FALSE)
    log_msg("   Mouse species: MouseRNAseqData (Benayoun et al. 2019 Genome Res)")
    log_msg("   Coverage: Immune, Neural, Endothelial, Mesenchymal, Hepatocyte, Cardiomyocyte and other major lineages (bulk RNA-seq deconvolution)")
  }
  log_msg(sprintf("   Number of reference broad labels: %d", length(unique(ref_broad$label.main))))
  log_msg(sprintf("   Label list: %s", paste(unique(ref_broad$label.main), collapse=", ")))

  log_msg("[L1-3] SingleR broad type annotation (per-cell)...")
  singler_broad <- SingleR(
    test            = sce_test,
    ref             = ref_broad,
    labels          = ref_broad$label.main,
    genes           = SINGLER_GENES,
    sd.thresh       = SINGLER_SD_THRESH,
    de.method       = SINGLER_DE_METHOD,
    quantile        = SINGLER_QUANTILE,
    fine.tune       = TRUE,
    aggr.ref        = TRUE,
    tune.thresh     = SINGLER_TUNE_THRESH,
    prune           = TRUE,
    assay.type.test = SINGLER_ASSAY,
    assay.type.ref  = SINGLER_ASSAY,
    BPPARAM         = BiocParallel::SerialParam()
  )
  pruned_flag <- pruneScores(singler_broad,
    nmads = PRUNE_N_MADS,
    min.diff.med  = PRUNE_MIN_DIFF_MED,
    min.diff.next = PRUNE_MIN_DIFF_NEXT)
  singler_broad$pruned.labels <- singler_broad$labels
  singler_broad$pruned.labels[pruned_flag] <- NA_character_
  broad_per_cell <- ifelse(is.na(singler_broad$pruned.labels),
                           "Unknown", singler_broad$pruned.labels)

  broad_label_cols <- colnames(singler_broad$scores)
  broad_score_per_cell <- sapply(seq_len(ncol(seurat_obj)), function(i) {
    lb <- singler_broad$pruned.labels[i]
    if (is.na(lb) || !lb %in% broad_label_cols) return(0)
    sc <- singler_broad$scores[i, lb]
    if (is.na(sc) || !is.finite(sc)) return(0)
    sc
  })

  get_mode  <- function(v) { v <- v[!is.na(v) & v != ""]; if (!length(v)) return("Unknown"); names(sort(table(v),decreasing=TRUE))[1] }
  get_med   <- function(v) { v <- v[is.finite(v) & v > 0]; if (!length(v)) return(0); median(v,na.rm=TRUE) }

  L1_raw_vec   <- sapply(cluster_levels, function(cl) get_mode(broad_per_cell[original_clusters==cl]))
  L1_score_vec <- sapply(cluster_levels, function(cl) get_med(broad_score_per_cell[original_clusters==cl]))
  L1_main_vec <- ifelse(
    L1_score_vec >= HPCA_CONF_THRESH & L1_raw_vec != "Unknown",
    L1_raw_vec,
    "Unknown"
  )

  log_msg(sprintf("\n  %s confidence threshold: %.2f (labels below this value are marked as Unknown)",
                  if (is_mouse) "MouseRNAseqData" else "HPCA", HPCA_CONF_THRESH))
  log_msg(sprintf("  %-8s %-30s %-10s → %-30s", "Cluster", "Reference dataset original labels", "Score", "Layer1 broad type"))
  for (i in seq_along(cluster_levels))
    log_msg(sprintf("  %-8s %-30s %-10.3f → %-30s",
                    cluster_levels[i], L1_raw_vec[i], L1_score_vec[i], L1_main_vec[i]))


  ref_name <- if (is_mouse) "MouseRNAseqData" else "HPCA"
  log_msg(sprintf("\n[L1-4] %s labels → marker file broad types fuzzy matching...", ref_name))
  hpca_unique_labels <- unique(L1_main_vec[L1_main_vec != "Unknown"])
  hpca_to_marker_map <- character(0)

  for (lb in hpca_unique_labels) {
    matched <- fuzzy_match_hpca(lb, marker_main_types, hardcode_map = ACTIVE_HARDCODE_MAP)
    hpca_to_marker_map[lb] <- if (is.na(matched)) lb else matched
    log_msg(sprintf("   %s '%s' → marker type '%s'%s",
                    ref_name, lb,
                    if (is.na(matched)) "(no match, using original label)" else matched,
                    if (is.na(matched)) " [!]" else ""))
  }
  hpca_to_marker_map["Unknown"] <- "Unknown"

  L1_mapped_vec <- sapply(L1_main_vec, function(lb) {
    if (lb %in% names(hpca_to_marker_map)) hpca_to_marker_map[lb] else lb
  })


  log_msg("\n========== [Layer 2] UCell → Subtype Annotation==========")

  ucell_expr <- GetAssayData(seurat_obj, assay=use_assay, layer="counts")


  L2_fine_vec  <- character(length(cluster_levels))
  L2_score_vec <- numeric(length(cluster_levels))

  if (length(fine_markers) > 0) {
    log_msg(sprintf("   Calculating UCell scores for %d fine types...", length(fine_markers)))
    ucell_fine_scores <- ScoreSignatures_UCell(
      matrix   = ucell_expr,
      features = fine_markers,
      maxRank  = UCELL_MAXRANK,
      BPPARAM  = BiocParallel::SerialParam()
    )


    fine_cluster_scores <- sapply(names(fine_markers), function(ft) {
      sapply(cluster_levels, function(cl) {
        col_nm <- paste0(ft, "_UCell")
        if (!col_nm %in% colnames(ucell_fine_scores)) return(0)
        sc <- ucell_fine_scores[original_clusters == cl, col_nm]
        sc <- sc[is.finite(sc)]
        if (!length(sc)) return(0)
        median(sc, na.rm=TRUE)
      })
    })


    for (i in seq_along(cluster_levels)) {
      cl        <- cluster_levels[i]
      mapped_main <- L1_mapped_vec[i]

      if (mapped_main == "Unknown") {
        L2_fine_vec[i]  <- "Unknown"
        L2_score_vec[i] <- 0
        log_msg(sprintf("  Cluster %2s: Layer1=Unknown → fine=Unknown", cl))
        next
      }


      candidates <- names(fine_to_main_db[fine_to_main_db == mapped_main])
      candidates <- candidates[candidates %in% colnames(fine_cluster_scores)]

      if (length(candidates) == 0) {

        L2_fine_vec[i]  <- mapped_main
        L2_score_vec[i] <- 0
        log_msg(sprintf("  Cluster %2s: %s 无亚型 marker → fine=%s",
                        cl, mapped_main, mapped_main))
        next
      }

      sub_vec  <- as.numeric(fine_cluster_scores[i, candidates, drop=TRUE])
      names(sub_vec) <- candidates
      sub_sorted <- sort(sub_vec, decreasing=TRUE)
      best_ft    <- names(sub_sorted)[1]
      best_sc    <- sub_sorted[1]
      second_sc  <- if (length(sub_sorted) >= 2) sub_sorted[2] else 0
      margin     <- best_sc - second_sc

      if (margin >= UCELL_MARGIN && best_sc > 0) {
        L2_fine_vec[i]  <- best_ft
        L2_score_vec[i] <- best_sc
        log_msg(sprintf("  Cluster %2s: [%s] → %s (UCell=%.3f, margin=%.3f)",
                        cl, mapped_main, best_ft, best_sc, margin))
      } else {
        L2_fine_vec[i]  <- mapped_main
        L2_score_vec[i] <- 0
        log_msg(sprintf("  Cluster %2s: [%s] 亚型竞争不足(margin=%.3f) → fine=%s",
                        cl, mapped_main, margin, mapped_main))
      }
    }
  } else {

    L2_fine_vec  <- L1_mapped_vec
    L2_score_vec <- rep(0, length(cluster_levels))
    log_msg("   No fine markers available, skipping Layer 2 UCell annotation.")
    }


  final_main_vec  <- L1_mapped_vec
  final_fine_vec  <- L2_fine_vec
  final_score_vec <- L1_score_vec

  log_msg("\n========== [Summary] Final Annotation ==========")
  log_msg(sprintf("  %-8s %-30s → %-38s  HPCA_score  UCell_fine_score",
                  "Cluster", "final_main(HPCA)", "final_fine(UCell)"))
  for (i in seq_along(cluster_levels))
    log_msg(sprintf("  %-8s %-30s → %-38s  %.3f       %.3f",
                    cluster_levels[i], final_main_vec[i], final_fine_vec[i],
                    L1_score_vec[i], L2_score_vec[i]))


  log_msg("\n[Steps] Mapping annotations to each cell...")
  cluster_to_main  <- setNames(final_main_vec,  cluster_levels)
  cluster_to_fine  <- setNames(final_fine_vec,  cluster_levels)
  cluster_to_score <- setNames(final_score_vec, cluster_levels)

  seurat_obj$celltype_main    <- unname(cluster_to_main[original_clusters])
  seurat_obj$celltype_fine    <- unname(cluster_to_fine[original_clusters])
  seurat_obj$annotation_score <- unname(cluster_to_score[original_clusters])

  stopifnot(!any(is.na(seurat_obj$celltype_fine)))
  stopifnot(!any(seurat_obj$celltype_fine == ""))
  log_msg("   Validation passed: celltype_main and celltype_fine are both non-empty")


  main_summary <- sort(table(seurat_obj$celltype_main), decreasing=TRUE)
  fine_summary <- sort(table(seurat_obj$celltype_fine), decreasing=TRUE)

  log_msg(paste0("\n", paste(rep("=",60), collapse="")))
  log_msg("Final Annotation Results")
  log_msg(paste(rep("=",60), collapse=""))
  log_msg(sprintf("Broad types %d, Subtypes %d",
                  length(main_summary), length(fine_summary)))

  log_msg("\nBroad type distribution")
  for (mt in names(main_summary))
    log_msg(sprintf("  %-30s %5d cells (%.1f%%)", mt, main_summary[mt],
                    100*main_summary[mt]/ncol(seurat_obj)))

  log_msg("\nSubtype distribution")
  for (ft in names(fine_summary))
    log_msg(sprintf("  %-42s %5d cells (%.1f%%)", ft, fine_summary[ft],
                    100*fine_summary[ft]/ncol(seurat_obj)))


  log_msg("\nGenerating visualization report...")
  main_levels <- sort(unique(seurat_obj$celltype_main))
  fine_levels <- sort(unique(seurat_obj$celltype_fine))
  main_cols <- setNames(safe_hue_pal(length(main_levels), h=c(15,375), l=65, c=100), main_levels)
  fine_cols <- setNames(safe_hue_pal(length(fine_levels), h=c(0,360),  l=70, c=90),  fine_levels)

  sapply(c(dirname(output_seurat), dirname(output_report), dirname(output_scores)),
         dir.create, showWarnings=FALSE, recursive=TRUE)

  pdf(output_report, width=16, height=8, onefile=TRUE)
  on.exit(try(dev.off(), silent=TRUE), add=TRUE)

  p1 <- DimPlot(seurat_obj, reduction="umap", group.by="seurat_clusters",
                label=TRUE, repel=TRUE, pt.size=0.5) +
    NoLegend() + labs(title="Clusters")
  p2 <- DimPlot(seurat_obj, reduction="umap", group.by="celltype_main",
                label=TRUE, repel=TRUE, cols=main_cols, pt.size=0.5) +
    NoLegend() + labs(title=sprintf("Major Cell Types [%s] (%s)",
                                    if (is_mouse) "MouseRNAseqData" else "HPCA",
                                    species_label))
  p3 <- DimPlot(seurat_obj, reduction="umap", group.by="celltype_fine",
                label=TRUE, repel=TRUE, cols=fine_cols, pt.size=0.5) +
    NoLegend() + labs(title=sprintf("Cell Subtypes [UCell] (%s)", species_label))
  print(p1 | p2 | p3)


  score_df <- data.frame(
    cluster   = factor(cluster_levels, levels=cluster_levels),
    hpca_score = L1_score_vec,
    fine_score = L2_score_vec,
    main_type  = final_main_vec,
    fine_type  = final_fine_vec,
    stringsAsFactors = FALSE
  )
  p_score <- ggplot(score_df, aes(x=cluster, y=hpca_score, fill=main_type)) +
    geom_col(width=0.7) +
    geom_hline(yintercept=HPCA_CONF_THRESH, linetype="dashed", color="red", linewidth=0.8) +
    annotate("text", x=1, y=HPCA_CONF_THRESH+0.01,
             label=sprintf("Confidence Threshold %.2f", HPCA_CONF_THRESH),
             hjust=0, size=3, color="red") +
    geom_text(aes(label=fine_type), angle=90, hjust=-0.05, size=2.2) +
    scale_fill_manual(values=main_cols) +
    labs(title=sprintf("%s Annotation Confidence Score per Cluster",
                       if (is_mouse) "MouseRNAseqData" else "HPCA"),
         subtitle=sprintf("SingleR Spearman r（%sReference set）; Red line = confidence threshold",
                          if (is_mouse) "MouseRNAseqData" else "HPCA"),
         x="Cluster", y="Reference Score", fill="Main Type") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    expand_limits(y=max(L1_score_vec, na.rm=TRUE)*1.4)
  print(p_score)


  main_counts <- as.data.frame(table(seurat_obj$celltype_main))
  colnames(main_counts) <- c("CellType","Count")
  main_counts$Percentage <- 100*main_counts$Count/sum(main_counts$Count)
  p_main <- ggplot(main_counts, aes(x=reorder(CellType,Count), y=Count, fill=CellType)) +
    geom_bar(stat="identity", show.legend=FALSE) +
    geom_text(aes(label=sprintf("%d (%.1f%%)", Count, Percentage)), hjust=-0.1, size=3) +
    coord_flip() + scale_fill_manual(values=main_cols) +
    labs(title=sprintf("Major Cell Types Distribution (%s)", species_label), x="", y="Cell Count") +
    theme_bw()
  print(p_main)


  fine_counts <- as.data.frame(table(seurat_obj$celltype_fine))
  colnames(fine_counts) <- c("CellType","Count")
  fine_counts$Percentage <- 100*fine_counts$Count/sum(fine_counts$Count)
  fine_counts <- fine_counts[order(-fine_counts$Count), ]
  if (nrow(fine_counts) > 30) fine_counts <- fine_counts[1:30, ]
  p_fine <- ggplot(fine_counts, aes(x=reorder(CellType,Count), y=Count, fill=CellType)) +
    geom_bar(stat="identity", show.legend=FALSE) +
    geom_text(aes(label=sprintf("%d (%.1f%%)", Count, Percentage)), hjust=-0.05, size=2.5) +
    coord_flip() + scale_fill_manual(values=fine_cols) +
    labs(title=sprintf("Cell Subtypes Distribution (%s)", species_label), x="", y="Cell Count") +
    theme_bw() + theme(axis.text.y=element_text(size=7))
  print(p_fine)


  match_df <- data.frame(
    cluster         = cluster_levels,
    hpca_raw_label  = L1_raw_vec,
    hpca_score      = round(L1_score_vec, 3),
    marker_main     = final_main_vec,
    final_fine      = final_fine_vec,
    ucell_fine_score = round(L2_score_vec, 3),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(match_df)))
    log_msg(sprintf("  Cluster %-4s HPCA='%-25s'(%.3f) → main='%-25s' fine='%-30s'(UCell=%.3f)",
                    match_df$cluster[i], match_df$hpca_raw_label[i], match_df$hpca_score[i],
                    match_df$marker_main[i], match_df$final_fine[i], match_df$ucell_fine_score[i]))

  if ("group" %in% colnames(seurat_obj@meta.data)) {
    groups <- unique(seurat_obj$group)
    groups <- groups[!is.na(groups)]
    for (g in groups) {
      sub_obj <- subset(seurat_obj, group == g)
      if (ncol(sub_obj) == 0) next
      pg1 <- DimPlot(sub_obj, reduction="umap", group.by="celltype_main",
                     cols=main_cols, label=TRUE, repel=TRUE, pt.size=0.5) +
        NoLegend() + labs(title=paste("Major Types -", g))
      pg2 <- DimPlot(sub_obj, reduction="umap", group.by="celltype_fine",
                     cols=fine_cols, label=TRUE, repel=TRUE, pt.size=0.5) +
        NoLegend() + labs(title=paste("Subtypes -", g))
      print(pg1 | pg2)
    }
  }

  dev.off()
  on.exit(NULL)


  fine_to_main <- setNames(final_main_vec, final_fine_vec)
  fine_to_main <- fine_to_main[!duplicated(names(fine_to_main))]

  annotation_results <- data.frame(
    cluster          = cluster_levels,
    main_type        = final_main_vec,
    fine_type        = final_fine_vec,
    hpca_raw_label   = L1_raw_vec,
    hpca_score       = L1_score_vec,
    ucell_fine_score = L2_score_vec,
    stringsAsFactors = FALSE
  )
  write.csv(annotation_results, output_scores, row.names=FALSE)

  seurat_obj@misc$annotation_info <- list(
    method         = "Hierarchical_Ref_UCell",
    architecture   = sprintf("Layer1: SingleR(%s) → main type; Layer2: UCell(marker内竞争) → fine type",
                             if (is_mouse) "MouseRNAseqData" else "HPCA"),
    reference_dataset = if (is_mouse) "MouseRNAseqData (Benayoun et al. 2019 Genome Res)"
                        else "HumanPrimaryCellAtlasData (Mabbott et al. 2013 BMC Genomics)",
    assay_used     = sprintf("%s (log-normalized for SingleR; raw counts for UCell)", use_assay),
    organism       = organism,
    hpca_conf_thresh = HPCA_CONF_THRESH,
    references = c(
      "Aran et al. (2019) Nature Immunology 20:163-172",
      "Mabbott et al. (2013) BMC Genomics 14:632",
      "Benayoun et al. (2019) Genome Res 29:1653-1665",
      "Andreatta & Carmona (2021) Comp Struct Biotechnol J 19:3796-3798",
      "Hao et al. (2021) Cell 184:3573-3587",
      "Luecken et al. (2022) Nature Methods 19:685-698",
      "Kikuta et al. (2025) Front Immunol 16:1614230"
    ),
    ucell_params   = list(maxRank=UCELL_MAXRANK, minMarkers=UCELL_MINMARKERS, margin=UCELL_MARGIN),
    singler_params = list(genes=SINGLER_GENES, de_method=SINGLER_DE_METHOD,
                          quantile=SINGLER_QUANTILE, prune_n_mads=PRUNE_N_MADS),
    n_main_types   = length(unique(seurat_obj$celltype_main)),
    n_fine_types   = length(unique(seurat_obj$celltype_fine)),
    date           = Sys.time()
  )
  seurat_obj@misc$fine_to_main_map   <- fine_to_main
  seurat_obj@misc$annotation_results <- annotation_results
  seurat_obj@misc$hpca_singler_result <- singler_broad
  seurat_obj$singler_labels          <- seurat_obj$celltype_fine

  saveRDS(seurat_obj, output_seurat, compress=FALSE)

  log_msg("\n=== Completion ===")
  log_msg(sprintf("assay %s | species %s | broad types %d | subtypes %d",
                  use_assay, species_label,
                  length(unique(seurat_obj$celltype_main)),
                  length(unique(seurat_obj$celltype_fine))))

}, error = function(e) {
  log_msg(paste("Error:", e$message), level="ERROR")
  traceback_info <- capture.output(traceback())
  for (line in traceback_info) log_msg(line, level="ERROR")
  stop(e)
})
