# rules/11_enrich.smk
rule run_pathway_enrichment:
    input:
        de_results = config["output"]["results"]["differential"] + "/spreadsheets/differential_expression_results.csv",
        final_annotated_seurat = config["output"]["results"]["annotation"] + "/seurat_final_annotated.rds"
    output:
        go_dir = directory(config["output"]["results"]["enrichment"] + "/GO"), 
        kegg_dir = directory(config["output"]["results"]["enrichment"] + "/KEGG"),
        enrichment_summary = config["output"]["results"]["enrichment"] + "/enrichment_summary.csv",
        enrichment_summary_txt = config["output"]["results"]["enrichment"] + "/enrichment_summary.txt",
        enrichment_rds = config["output"]["results"]["enrichment"] + "/enrichment_results_all.rds",
        unmapped_genes = config["output"]["results"]["enrichment"] + "/unmapped_genes.csv",
        enrichment_log = config["output"]["results"]["enrichment"] + "/pathway_enrichment.log"
    params:
        organism = config["organism"],
        mouse_db_path = config["downstream"]["pathway_enrichment"]["databases"]["mouse_db_path"],
        human_db_path = config["downstream"]["pathway_enrichment"]["databases"]["human_db_path"],
        pval_cutoff = config["downstream"]["pathway_enrichment"]["parameters"]["pval_cutoff"],
        qval_cutoff = config["downstream"]["pathway_enrichment"]["parameters"]["qval_cutoff"],
        logfc_threshold = config["downstream"]["pathway_enrichment"]["parameters"]["logfc_threshold"],
        min_gene_size = config["downstream"]["pathway_enrichment"]["parameters"]["min_gene_size"],
        max_gene_size = config["downstream"]["pathway_enrichment"]["parameters"]["max_gene_size"]
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["enrichment"] + "/logs/run_pathway_enrichment.log"
    threads:
        config["resources"]["downstream"]["pathway_enrichment"]["threads"]
    resources:
        mem_gb = config["resources"]["downstream"]["pathway_enrichment"]["mem_gb"]
    script:
        config["script_paths"]["enrichment_scripts"] + "/run_pathway_enrichment.R"
rule visualize_enrichment:
    input:
        enrichment_rds = rules.run_pathway_enrichment.output.enrichment_rds,
        enrichment_summary = rules.run_pathway_enrichment.output.enrichment_summary
    output:
        go_viz_dir = directory(config["output"]["results"]["enrichment"] + "/GO_visualization"),
        kegg_viz_dir = directory(config["output"]["results"]["enrichment"] + "/KEGG_visualization"),
        go_done = config["output"]["results"]["enrichment"] + "/GO_visualization.done",
        kegg_done = config["output"]["results"]["enrichment"] + "/KEGG_visualization.done"
    params:
        top_n_terms = config["downstream"]["pathway_enrichment"]["visualization"]["top_n_terms"],
        plot_width = config["downstream"]["pathway_enrichment"]["visualization"]["plot_width"],
        plot_height = config["downstream"]["pathway_enrichment"]["visualization"]["plot_height"]
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["enrichment"] + "/logs/visualize_enrichment.log"
    threads:
        config["resources"]["downstream"]["pathway_enrichment"]["threads"]
    resources:
        mem_gb = config["resources"]["downstream"]["pathway_enrichment"]["mem_gb"]
    script:
        config["script_paths"]["enrichment_scripts"] + "/visualize_enrichment.R"
rule enrichment_all:
    input:
        rules.run_pathway_enrichment.output.enrichment_summary,
        rules.run_pathway_enrichment.output.enrichment_rds,
        rules.visualize_enrichment.output.go_viz_dir,
        rules.visualize_enrichment.output.kegg_viz_dir

        