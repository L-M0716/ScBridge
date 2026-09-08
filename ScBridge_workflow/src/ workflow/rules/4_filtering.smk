# rules/4_filtering.smk
with open(config["input"]["samples_table"]) as f:
    next(f)
    SAMPLES = [line.strip().split("\t")[0] for line in f if line.strip()]
rule filter_cells_one_sample:
    input:
        seurat_with_doublets = config["output"]["results"]["doublets"] + "/{sample_id}/seurat_with_doublets.rds"
    output:
        filtered_seurat = config["output"]["results"]["filtering"] + "/{sample_id}/seurat_filtered.rds",
        filter_stats = config["output"]["results"]["filtering"] + "/{sample_id}/filtering_stats.csv",
        filter_report = config["output"]["results"]["filtering"] + "/{sample_id}/filtering_report.pdf"
    params:
        min_genes = config["qc"]["min_genes"],
        max_genes = config["qc"]["max_genes"],
        min_counts = config["qc"]["min_counts"],
        max_mito = config["qc"]["max_mito"],
        min_cells_per_gene = config["qc"]["gene_filtering"]["min_cells"],
        exclude_genes = ",".join(config["qc"]["gene_filtering"]["exclude_genes"])
    conda:
        config["conda"]["sctest"]
    log:
        config["output"]["results"]["filtering"] + "/logs/{sample_id}_filtering.log"
    threads:
        config["resources"]["filtering"]["threads"]
    resources:
        mem_gb = config["resources"]["filtering"]["mem_gb"]
    script:
        config["script_paths"]["filtering_scripts"] + "/filter_cells_one_sample.R"

rule filtering_all:
    input:
        expand(rules.filter_cells_one_sample.output.filtered_seurat, sample_id=SAMPLES)