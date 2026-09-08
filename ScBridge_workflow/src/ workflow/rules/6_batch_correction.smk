# rules/6_batch_correction.smk
rule run_batch_correction:
    input:
        seurat_normalized = config["output"]["results"]["normalization"] + "/seurat_merged_normalized.rds"
    output:
        corrected = config["output"]["results"]["batch_correction"] + "/batch_corrected.rds",
        umap_plot = config["output"]["results"]["batch_correction"] + "/batch_correction_umap.pdf",
        pca_plot = config["output"]["results"]["batch_correction"] + "/batch_correction_pca.pdf"
    params:
        method = config["integration"]["method"],
        batch_key = config["integration"]["batch_key"],
        seed = config["normalization"]["seed"],
        n_pcs = config["integration"]["n_pcs"]
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["batch_correction"] + "/logs/batch_correction.log"
    threads:
        config["resources"]["batch_correction"]["correction"]["threads"]
    resources:
        mem_gb = config["resources"]["batch_correction"]["correction"]["mem_gb"]
    script:
        config["script_paths"]["batch_correction_scripts"] + "/correct_batches.R" 
rule batch_correction_all:
    input:
        rules.run_batch_correction.output.corrected,
        rules.run_batch_correction.output.umap_plot