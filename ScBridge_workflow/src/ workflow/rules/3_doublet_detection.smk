# rules/3_doublet_detection.smk

with open(config["input"]["samples_table"]) as f:
    next(f)
    SAMPLES = [line.strip().split("\t")[0] for line in f if line.strip()]

rule detect_doublets_one_sample:
    input:
        seurat_qc = config["output"]["results"]["qc"] + "/{sample_id}/seurat_qc.rds"
    output:
        seurat_with_doublets = config["output"]["results"]["doublets"] + "/{sample_id}/seurat_with_doublets.rds",
        doublet_stats = config["output"]["results"]["doublets"] + "/{sample_id}/doublet_stats.csv"
    params:
        expected_rate = config["qc"]["expected_doublet_rate"],
        method = config["qc"]["doublet_method"],
        seed = config["normalization"]["seed"]
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["doublets"] + "/logs/{sample_id}_detection.log"
    threads:
        config["resources"]["doublet_detection"]["threads"]
    resources:
        mem_gb = config["resources"]["doublet_detection"]["mem_gb"]
    script:
        config["script_paths"]["doublet_scripts"] + "/detect_doublets.R"

rule doublet_detection_all:
    input:
        expand(rules.detect_doublets_one_sample.output.seurat_with_doublets, sample_id=SAMPLES)