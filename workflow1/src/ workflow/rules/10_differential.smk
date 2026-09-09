# rules/10_differential.smk
rule run_differential_expression:
    input:
        annotated_seurat = config["output"]["results"]["annotation"] + "/seurat_final_annotated.rds",
        sample_table = config["input"]["samples_table"]
    output:
        de_results_table = config["output"]["results"]["differential"] + "/spreadsheets/differential_expression_results.csv"
    params:
        group_by_column = config["differential"]["group_by_column"],
        control_group = config["differential"]["control_group"],
        logfc_threshold = config["differential"]["logfc_threshold"],
        min_cells_group = config["differential"]["min_cells_group"]
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["differential"] + "/logs/run_differential_expression.log"
    threads:
        config["resources"]["differential"]["run"]["threads"]
    resources:
        mem_gb = config["resources"]["differential"]["run"]["mem_gb"]
    script:
        config["script_paths"]["differential_scripts"] + "/run_differential_expression.R"
rule visualize_differential_expression:
    input:
        de_results_table = rules.run_differential_expression.output.de_results_table,
        annotated_seurat = rules.run_differential_expression.input.annotated_seurat
    output:
        volcano_plots = config["output"]["results"]["differential"] + "/reports/de_volcano_plots.pdf",
        de_heatmap = config["output"]["results"]["differential"] + "/reports/de_heatmap.pdf",
        top_de_genes_table = config["output"]["results"]["differential"] + "/spreadsheets/top_de_genes_summary.csv"
    params:
        top_n = config["differential"]["visualization"]["top_n"],
        adj_pval_cutoff = config["differential"]["visualization"]["adj_pval_cutoff"],
        vis_logfc_threshold = config["differential"]["visualization"]["logfc_threshold"]
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["differential"] + "/logs/visualize_differential_expression.log"
    threads:
        config["resources"]["differential"]["viz"]["threads"]
    resources:
        mem_gb = config["resources"]["differential"]["viz"]["mem_gb"]
    script:
        config["script_paths"]["differential_scripts"] + "/visualize_differential_expression.R"
rule differential_all:
    input:
        rules.run_differential_expression.output.de_results_table,
        rules.visualize_differential_expression.output.volcano_plots,
        rules.visualize_differential_expression.output.de_heatmap,
        rules.visualize_differential_expression.output.top_de_genes_table