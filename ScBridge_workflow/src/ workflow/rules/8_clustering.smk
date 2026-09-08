# rules/8_clustering.smk
rule run_clustering:
    input:
        purified_seurat = config["output"]["results"]["feature_selection"] + "/seurat_feature.rds",
        sample_table = config["input"]["samples_table"]
    output:
        clustered_seurat = config["output"]["results"]["clustered"] + "/seurat_clustered.rds",
        umap_plot_preview = config["output"]["results"]["clustered"] + "/clustering_preview_umap.pdf"
    params:
        n_pcs = config["integration"]["n_pcs"],
        reduction_method = config["integration"]["method"],
        n_neighbors = config["clustering"]["n_neighbors"],
        resolution = config["clustering"]["resolution"]["default"],
        algorithm = config["clustering"]["algorithm"],
        min_dist = config["clustering"]["umap"]["min_dist"],
        spread = config["clustering"]["umap"]["spread"],
        seed = config["normalization"]["seed"]
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["clustered"] + "/logs/clustering.log"
    threads:
        config["resources"]["clustering"]["pca"]["threads"]
    resources:
        mem_gb = config["resources"]["clustering"]["pca"]["mem_gb"]
    script:
        config["script_paths"]["clustering_scripts"] + "/run_clusters.R"

### 元规则: 
rule clustering_all:
    input:
        rules.run_clustering.output.clustered_seurat,
        rules.run_clustering.output.umap_plot_preview