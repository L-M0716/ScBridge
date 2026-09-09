# scripts/normalization_scripts/normalize_one_sample.R
suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
  library(glmGamPoi)
  library(future)
})

input_file      <- snakemake@input[["seurat_filtered"]]
cell_cycle_genes_path <- snakemake@input[["cell_cycle_genes"]]
output_file     <- snakemake@output[["seurat_normalized"]]
confounders     <- snakemake@params[["vars_to_regress"]]
min_cells_per_gene <- snakemake@params[["min_cells_per_gene"]]
return_only_var_genes <- as.logical(snakemake@params[["return_only_var_genes"]])
seed            <- snakemake@params[["seed"]]
log_file        <- snakemake@log[[1]]
sample_id       <- snakemake@wildcards[["sample_id"]]

dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_msg <- function(msg, level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] [%s] [%s] %s\n", ts, level, sample_id, msg), file = log_file, append = TRUE)
}
tryCatch({
  set.seed(seed)
  plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  log_msg("Script started (single sample normalization).")
  obj <- readRDS(input_file)
  log_msg(sprintf("Loaded sample %s, cells=%d, genes=%d", sample_id, ncol(obj), nrow(obj)))
  
  DefaultAssay(obj) <- "RNA"
  obj[["RNA"]] <- tryCatch(JoinLayers(obj[["RNA"]]), error = function(e) obj[["RNA"]])
  log_msg("NormalizeData -> Generate RNA data layer")
  obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)
  
  if (!file.exists(cell_cycle_genes_path)) {
    log_msg("Cell cycle gene file does not exist, skipping scoring", "WARNING")
  } else {
    cc <- tryCatch(readRDS(cell_cycle_genes_path), error = function(e) NULL)
    if (is.null(cc) || !all(c("s.genes","g2m.genes") %in% names(cc))) {
      log_msg("Cell cycle gene file has incorrect format, skipping scoring", "WARNING")
    } else {
      s_matched  <- intersect(cc$s.genes, rownames(obj))
      g2m_matched<- intersect(cc$g2m.genes, rownames(obj))
      if (length(s_matched) >= 5 && length(g2m_matched) >= 5) {
        log_msg(sprintf("CellCycleScoring: S=%d, G2M=%d", length(s_matched), length(g2m_matched)))
        obj <- CellCycleScoring(obj, s.features = s_matched, g2m.features = g2m_matched,
                                assay = "RNA", set.ident = FALSE, verbose = FALSE)
      } else {
        log_msg(sprintf("Cell cycle genes insufficient (S=%d, G2M=%d), skipping scoring", length(s_matched), length(g2m_matched)), "WARNING")
      }
    }
  }
  
  log_msg(sprintf("SCTransform（vars_to_regress=%s）", ifelse(length(confounders)==0,"无", paste(confounders,collapse=","))))
  obj <- SCTransform(
    obj,
    assay = "RNA",
    new.assay.name = "SCT",
    vars.to.regress = confounders,
    method = "glmGamPoi",
    min_cells = min_cells_per_gene,
    variable.features.n = 3000,
    return.only.var.genes = return_only_var_genes,
    vst.flavor = "v2",
    seed.use = seed,
    verbose = FALSE
  )
  DefaultAssay(obj) <- "SCT"
  log_msg(sprintf("SCT completed, HVG count: %d", length(VariableFeatures(obj))))
  
  saveRDS(obj, file = output_file)
  log_msg(paste("Normalized object saved to:", output_file))
  log_msg("Script executed successfully.")
}, error = function(e) {
  log_msg(paste("ERROR:", e$message), level = "ERROR")
  stop(e)
})
