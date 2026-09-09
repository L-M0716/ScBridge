# rules/13_cellchat.smk
rule cellchat_inference:
    input:
        annotated = config["output"]["results"]["annotation"] + "/seurat_final_annotated.rds"
    output:
        cellchat_obj = config["output"]["results"]["cellchat"] + "/cellchat_object.rds",
        net_df = config["output"]["results"]["cellchat"] + "/spreadsheets/cell_communication_network.csv",
        lr_pairs = config["output"]["results"]["cellchat"] + "/spreadsheets/LR_pairs_all.csv"
    params:
        db_type = config["downstream"]["cellchat"]["db"],
        min_cells = config["downstream"]["cellchat"]["min_cells"],
        cell_type_column = config["downstream"]["cellchat"]["cell_type_column"],
        organism = config["organism"]
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["cellchat"] + "/logs/cellchat_inference.log"
    threads:
        config["resources"]["downstream"]["cellchat"]["threads"]
    resources:
        mem_mb = config["resources"]["downstream"]["cellchat"]["mem_gb"]*1024
    script:
        config["script_paths"]["cellchat_scripts"] + "/run_cellchat_inference.R"
rule cellchat_visualization_network:
    input:
        cellchat_obj = rules.cellchat_inference.output.cellchat_obj
    output:
        plot_dir = directory(config["output"]["results"]["cellchat"] + "/plots/network"),
        circle_plot = config["output"]["results"]["cellchat"] + "/plots/network/interaction_circle.pdf",
        heatmap_count = config["output"]["results"]["cellchat"] + "/plots/network/interaction_heatmap_count.pdf",
        heatmap_weight = config["output"]["results"]["cellchat"] + "/plots/network/interaction_heatmap_weight.pdf"
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["cellchat"] + "/logs/cellchat_viz_network.log"
    threads:
        config["resources"]["downstream"]["cellchat"]["threads"]
    resources:
        mem_mb = config["resources"]["downstream"]["cellchat"]["mem_gb"]*1024
    script:
        config["script_paths"]["cellchat_scripts"] + "/run_cellchat_viz_network.R"


rule cellchat_visualization_pathway:
    input:
        cellchat_obj = rules.cellchat_inference.output.cellchat_obj
    output:
        plot_dir = directory(config["output"]["results"]["cellchat"] + "/plots/pathways"),
        pathway_summary = config["output"]["results"]["cellchat"] + "/spreadsheets/pathway_summary.csv",
        signaling_role = config["output"]["results"]["cellchat"] + "/plots/pathways/signaling_role_scatter.pdf"
    params:
        pathways_of_interest = config["downstream"]["cellchat"]["pathways_of_interest"],
        top_n_pathways = config["downstream"]["cellchat"]["top_n_pathways"]
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["cellchat"] + "/logs/cellchat_viz_pathway.log"
    threads:
        config["resources"]["downstream"]["cellchat"]["threads"]
    resources:
        mem_mb = config["resources"]["downstream"]["cellchat"]["mem_gb"]*1024
    script:
        config["script_paths"]["cellchat_scripts"] + "/run_cellchat_viz_pathway.R"


rule cellchat_pattern_analysis:
    input:
        cellchat_obj = rules.cellchat_inference.output.cellchat_obj
    output:
        cellchat_final = config["output"]["results"]["cellchat"] + "/cellchat_final.rds",
        outgoing_pattern = config["output"]["results"]["cellchat"] + "/plots/patterns/outgoing_pattern.pdf",
        incoming_pattern = config["output"]["results"]["cellchat"] + "/plots/patterns/incoming_pattern.pdf",
        signaling_role_heatmap = config["output"]["results"]["cellchat"] + "/plots/patterns/signaling_role_heatmap.pdf",
        centrality_scores = config["output"]["results"]["cellchat"] + "/spreadsheets/centrality_scores.csv"
    params:
        n_patterns_outgoing = config["downstream"]["cellchat"]["n_patterns_outgoing"],
        n_patterns_incoming = config["downstream"]["cellchat"]["n_patterns_incoming"]
    conda:
        config["conda"]["doublet_detect"]
    log:
        config["output"]["results"]["cellchat"] + "/logs/cellchat_pattern.log"
    threads:
        config["resources"]["downstream"]["cellchat"]["threads"]
    resources:
        mem_mb = config["resources"]["downstream"]["cellchat"]["mem_gb"]*1024
    script:
        config["script_paths"]["cellchat_scripts"] + "/run_cellchat_pattern.R"
rule cellchat_all:
    input:
        rules.cellchat_inference.output.cellchat_obj,
        rules.cellchat_visualization_network.output.circle_plot,
        rules.cellchat_visualization_pathway.output.signaling_role,
        rules.cellchat_pattern_analysis.output.cellchat_final
    default_target: True