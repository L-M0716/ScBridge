# rules/3.5_doublet_visualization.smk

with open(config["input"]["samples_table"]) as f:
    next(f)
    SAMPLES = [line.strip().split("\t")[0] for line in f if line.strip()]

rule visualize_doublets_one_sample:
    input:
        seurat_with_doublets = config["output"]["results"]["doublets"] + "/{sample_id}/seurat_with_doublets.rds"
    output:
        umap_plot_single = temp(config["output"]["results"]["doublets"] + "/{sample_id}_doublet_umap.pdf")
    params:
        confounders = config["normalization"]["vars_to_regress"],
        n_pcs = 15
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["doublets"] + "/logs/{sample_id}_visualization.log"
    threads:
        config["resources"]["normalization"]["sctransform"]["threads"]
    resources:
        mem_gb = config["resources"]["normalization"]["sctransform"]["mem_gb"]
    script:
        config["script_paths"]["doublet_scripts"] + "/visualize_doublets.R"

rule aggregate_doublet_plots:
    input:
        single_plots = expand(rules.visualize_doublets_one_sample.output.umap_plot_single, sample_id=SAMPLES)
    output:
        aggregated_plot = config["output"]["results"]["doublets"] + "/all_samples_doublet_umap.pdf"
    log:
        config["output"]["results"]["doublets"] + "/logs/aggregate_visualization.log"
    conda:
        config["conda"]["doublet_detect"]
    shell:
        "gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile={output.aggregated_plot} {input.single_plots}"

rule doublet_visualization_all:
    input:
        rules.aggregate_doublet_plots.output.aggregated_plot