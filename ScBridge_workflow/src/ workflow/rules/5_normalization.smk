# rules/5_normalization.smk

with open(config["input"]["samples_table"]) as f:
    next(f)
    SAMPLES = [line.strip().split("\t")[0] for line in f if line.strip()]


def filtered_rds_input(wildcards):
    filtered_by_sample = config.get("input", {}).get("filtered_rds_by_sample", {})
    if filtered_by_sample:
        sample_id = wildcards.sample_id
        if sample_id not in filtered_by_sample:
            raise ValueError(
                f"No filtered RDS path found for sample '{sample_id}' "
                "in config['input']['filtered_rds_by_sample']"
            )
        return filtered_by_sample[sample_id]

    return config["output"]["results"]["filtering"] + f"/{wildcards.sample_id}/seurat_filtered.rds"


rule normalize_one_sample:
    input:
        seurat_filtered = filtered_rds_input,
        cell_cycle_genes = config["input"]["cell_cycle_genes"]
    output:
        seurat_normalized = temp(config["output"]["results"]["normalization"] + "/{sample_id}/seurat_normalized.rds")
    params:
        vars_to_regress = config["normalization"]["vars_to_regress"],
        return_only_var_genes = False,
        seed = config["normalization"]["seed"],
        min_cells_per_gene = config["qc"]["gene_filtering"]["min_cells"]
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["normalization"] + "/logs/{sample_id}_normalization.log"
    threads:
        config["resources"]["normalization"]["sctransform"]["threads"]
    resources:
        mem_gb = config["resources"]["normalization"]["sctransform"]["mem_gb"]
    script:
        config["script_paths"]["normalization_scripts"] + "/run_sctransform.R"


rule aggregate_normalized_data:
    input:
        normalized_list = expand(rules.normalize_one_sample.output.seurat_normalized, sample_id=SAMPLES)
    output:
        merged_normalized_seurat = config["output"]["results"]["normalization"] + "/seurat_merged_normalized.rds",
        qc_figures = config["output"]["results"]["normalization"] + "/normalization_qc_figures.pdf",
        hvg_list = config["output"]["results"]["normalization"] + "/consensus_hvg_list.txt"
    params:
        n_variable_genes = config["normalization"]["n_variable_genes"]
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["normalization"] + "/logs/aggregate_normalization.log"
    threads:
        config["resources"]["normalization"]["sctransform"]["threads"]
    resources:
        mem_gb = config["resources"]["normalization"]["sctransform"]["mem_gb"]
    script:
        config["script_paths"]["normalization_scripts"] + "/aggregate_normalized.R"


rule normalization_all:
    input:
        rules.aggregate_normalized_data.output.merged_normalized_seurat
