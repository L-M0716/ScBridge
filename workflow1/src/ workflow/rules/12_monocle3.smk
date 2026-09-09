# rules/12_monocle3.smk
rule run_monocle3_analysis:
    input:
        annotated = config["output"]["results"]["annotation"] + "/seurat_final_annotated.rds"
    output:
        cds_object = config["output"]["results"]["monocle3"] + "/cds_object.rds",
        trajectory_summary = config["output"]["results"]["monocle3"] + "/trajectory_summary.rds",
        trajectory_summary_txt = config["output"]["results"]["monocle3"] + "/trajectory_summary.txt",
        trajectory_genes = config["output"]["results"]["monocle3"] + "/trajectory_dependent_genes.csv",
        gene_modules = config["output"]["results"]["monocle3"] + "/gene_modules.csv",
        gene_modules_summary = config["output"]["results"]["monocle3"] + "/gene_modules_summary.csv",
        plot_dir = directory(config["output"]["results"]["monocle3"] + "/plot")
    params:
        clusters_to_use = config["downstream"]["trajectory"]["clusters_to_use"],
        root_cell_type = config["downstream"]["trajectory"]["root_cell_type"],
        root_cluster = config["downstream"]["trajectory"]["root_cluster"],
        num_dim = config["downstream"]["trajectory"]["num_dim"],
        resolution = config["downstream"]["trajectory"]["resolution"],
        k = config["downstream"]["trajectory"]["k"],
        q_value_threshold = config["downstream"]["trajectory"]["q_value_threshold"]
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["monocle3"] + "/logs/monocle3.log"
    threads: 
        config["resources"]["downstream"]["monocle3"]["threads"]
    resources:
        mem_gb = config["resources"]["downstream"]["monocle3"]["mem_gb"]
    script:
        config["script_paths"]["monocle3_script"] + "/monocle3_analysis.R"
rule monocle3_all:
    input:
        rules.run_monocle3_analysis.output.cds_object,
        rules.run_monocle3_analysis.output.trajectory_summary,
        rules.run_monocle3_analysis.output.trajectory_summary_txt,
        rules.run_monocle3_analysis.output.trajectory_genes,
        rules.run_monocle3_analysis.output.gene_modules,
        rules.run_monocle3_analysis.output.gene_modules_summary,
        rules.run_monocle3_analysis.output.plot_dir
