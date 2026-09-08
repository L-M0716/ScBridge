# rules/7_feature_selection.smk
rule visualize_and_purify_features:
    input:
        corrected_seurat = config["output"]["results"]["batch_correction"] + "/batch_corrected.rds"
    output:
        feature_seurat = config["output"]["results"]["feature_selection"] + "/seurat_feature.rds",
        feature_plot = config["output"]["results"]["feature_selection"] + "/hvg_visualization_report.pdf"
    log:
        config["output"]["results"]["feature_selection"] + "/logs/visualize_features.log"
    threads:
        config["resources"]["default_threads"]
    conda:
        config["conda"]["doublet_detect"]
    resources:
        mem_gb = config["resources"]["default_mem_gb"]
    script:
        config["script_paths"]["feature_selection_scripts"] + "/feature_selection_visual.R"
rule feature_selection_all:
    input:
        rules.visualize_and_purify_features.output.feature_seurat,
        rules.visualize_and_purify_features.output.feature_plot