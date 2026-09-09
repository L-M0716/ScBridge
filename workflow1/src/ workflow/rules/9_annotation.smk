# rules/9_annotation.smk
rule auto_annotation:
    input:
        clustered_seurat = config["output"]["results"]["clustered"] + "/seurat_clustered.rds",
        markerlist = config["database_resources"]["reference_databases"]["markerlist"]
    output:
        auto_annotated = config["output"]["results"]["annotation"] + "/seurat_auto_annotated.rds",
        annotation_report = config["output"]["results"]["annotation"] + "/reports/auto_annotation_report.pdf",
        confidence_scores = config["output"]["results"]["annotation"] + "/spreadsheets/auto_annotation_scores.csv"
    params:
        organism = config["organism"]
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["annotation"] + "/logs/auto_annotation.log"
    threads:
        config["resources"]["annotation"]["auto"]["threads"]
    resources:
        mem_gb = config["resources"]["annotation"]["auto"]["mem_gb"]
    script:
        config["script_paths"]["annotation_scripts"] + "/run_auto_annotation.R"
rule analyze_and_visualize_annotation:
    input:
        auto_annotated = rules.auto_annotation.output.auto_annotated
    output:
        final_annotated_seurat = config["output"]["results"]["annotation"] + "/seurat_final_annotated.rds",
        umap_annotated = config["output"]["results"]["annotation"] + "/reports/umap_final_annotated.pdf",
        composition_plot = config["output"]["results"]["annotation"] + "/reports/celltype_composition.pdf",
        composition_table = config["output"]["results"]["annotation"] + "/spreadsheets/celltype_composition.csv",
        marker_heatmap = config["output"]["results"]["annotation"] + "/reports/marker_gene_heatmap.pdf"
    params:
        logfc_threshold = config["differential"]["logfc_threshold"],
        min_pct = config["differential"]["min_pct"]
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["annotation"] + "/logs/analyze_annotation.log"
    threads:
        config["resources"]["annotation"]["auto"]["threads"]
    resources:
        mem_gb = config["resources"]["annotation"]["auto"]["mem_gb"]
    script:
        config["script_paths"]["annotation_scripts"] + "/analyze_annotation.R"
rule annotation_all:
    input:
        rules.analyze_and_visualize_annotation.output.final_annotated_seurat,
        rules.analyze_and_visualize_annotation.output.umap_annotated
