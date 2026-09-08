# rules/2_quality_control.smk
with open(config["input"]["samples_table"]) as f:
    next(f)
    SAMPLES = [line.strip().split("\t")[0] for line in f if line.strip()]
rule quality_control_all:
    input:
        config["output"]["results"]["qc"] + "/aggregated_qc_report.pdf"
 
matrix_input_root = (
    config["input"]["matrix_dir"]
    if config.get("runtime", {}).get("input_mode") == "matrix"
    else config["output"]["results"]["base"] + "/read_data"
)

rule load_single_cell_data:
    input:
        matrix_dir = matrix_input_root + "/{sample_id}/filtered_feature_bc_matrix",
        sample_table = config["input"]["samples_table"]
    output:
        seurat_obj = config["output"]["results"]["qc"] + "/{sample_id}/seurat_raw.rds"
    params:
        mito_pattern = config["qc"]["mito_pattern"]
    conda:
        config["conda"]["sctest"]
    log:
        config["output"]["results"]["qc"] + "/logs/{sample_id}_load_data.log"
    threads:
        config["resources"]["default_threads"]
    resources:
        mem_gb = config["resources"]["default_mem_gb"]
    script:
        config["script_paths"]["io_scripts"] + "/load_10x_data.R"

rule calculate_qc_metrics:
    input:
        seurat_raw = rules.load_single_cell_data.output.seurat_obj
    output:
        qc_table = config["output"]["results"]["qc"] + "/{sample_id}/qc_metrics.csv",
        seurat_qc = config["output"]["results"]["qc"] + "/{sample_id}/seurat_qc.rds"
    params:
        min_genes = config["qc"]["min_genes"],
        max_genes = config["qc"]["max_genes"],
        min_counts = config["qc"]["min_counts"],
        max_mito = config["qc"]["max_mito"]
    conda:
        config["conda"]["sctest"]
    log:
        config["output"]["results"]["qc"] + "/logs/{sample_id}_calculate_qc.log"
    threads:
        config["resources"]["default_threads"]
    resources:
        mem_gb = config["resources"]["default_mem_gb"]
    script:
        config["script_paths"]["qc_scripts"] + "/calculate_qc_metrics.R"

rule aggregate_qc_plots:
    input:
        seurat_qc_list = expand(rules.calculate_qc_metrics.output.seurat_qc, sample_id=SAMPLES)
    output:
        aggregated_report = config["output"]["results"]["qc"] + "/aggregated_qc_report.pdf"
    conda:
        config["conda"]["sctest"]
    log:
        config["output"]["results"]["qc"] + "/logs/aggregate_qc_plots.log"
    threads:
        config["resources"]["default_threads"]
    resources:
        mem_gb = config["resources"]["default_mem_gb"]
    script:
        config["script_paths"]["qc_scripts"] + "/aggregate_qc_plots.R"