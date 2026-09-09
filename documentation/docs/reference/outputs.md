<div class="reference-tables" markdown="1">

# Output Reference

The tables below list analysis outputs relative to `RESULTS_DIR`, specified with `-R`.
The requested target and its dependencies determine which outputs are generated.

## Core results

| Stage | Files | Contents |
|---|---|---|
| Read QC | `read_data/runs/<run_id>/` | FastQC reports and read-processing files |
| Cell QC | <code>qc/<wbr>&lt;sample_<wbr>id&gt;/<wbr>qc_<wbr>metrics.<wbr>csv</code> | Per-cell QC metrics and status |
| QC report | <code>qc/<wbr>aggregated_<wbr>qc_<wbr>report.<wbr>pdf</code> | Cross-sample QC plots |
| Doublets | <code>doublets/<wbr>&lt;sample_<wbr>id&gt;/<wbr>doublet_<wbr>stats.<wbr>csv</code> | Sample doublet summary |
| Filtering | <code>filtering/<wbr>&lt;sample_<wbr>id&gt;/<wbr>seurat_<wbr>filtered.<wbr>rds</code> | Filtered Seurat object |
| Normalization | <code>normalization/<wbr>seurat_<wbr>merged_<wbr>normalized.<wbr>rds</code> | Merged normalized object |
| Integration | <code>batch_<wbr>correction/<wbr>batch_<wbr>corrected.<wbr>rds</code> | Integrated/corrected object |
| Features | <code>feature_<wbr>selection/<wbr>seurat_<wbr>feature.<wbr>rds</code> | Object passed to clustering |
| Clustering | <code>clustered/<wbr>seurat_<wbr>clustered.<wbr>rds</code> | Clustered object with UMAP |
| Annotation | <code>annotation/<wbr>seurat_<wbr>final_<wbr>annotated.<wbr>rds</code> | Final broad and fine cell labels |

## Downstream results

| Analysis | Files | Contents |
|---|---|---|
| Differential expression | <code>differential/<wbr>spreadsheets/<wbr>differential_<wbr>expression_<wbr>results.<wbr>csv</code> | Pseudobulk statistics |
| Differential reports | `differential/reports/` | Volcano plots and heatmap |
| Enrichment | <code>enrichment/<wbr>enrichment_<wbr>results_<wbr>all.<wbr>rds</code> | Collected enrichment objects |
| Enrichment summary | <code>enrichment/<wbr>enrichment_<wbr>summary.<wbr>csv</code> | Results and skip statuses |
| Enrichment plots | <code>enrichment/<wbr>GO_<wbr>visualization/<wbr></code>, <code>enrichment/<wbr>KEGG_<wbr>visualization/<wbr></code> | GO/KEGG figures |
| Trajectory | `monocle3/cds_object.rds` | Monocle3 object |
| Trajectory genes | <code>monocle3/<wbr>trajectory_<wbr>dependent_<wbr>genes.<wbr>csv</code> | Graph-associated genes retained below the q-value threshold |
| Trajectory plots | `monocle3/plot/` | Trajectory and pseudotime views |
| Communication | <code>cellchat/<wbr>cellchat_<wbr>object.<wbr>rds</code> | CellChat object or status object |
| Communication tables | `cellchat/spreadsheets/` | Network, ligand-receptor, pathway, and centrality tables |
| Communication plots | `cellchat/plots/` | Network, pathway, and pattern figures |

## Logs and runtime records

- `reports/standardization_report.txt` records input standardization.
- `runtime/` contains generated configuration and input manifests.
- Stage-specific `logs/` directories contain diagnostic messages.
- Per-sample normalized RDS files are declared as temporary outputs. Snakemake removes them after all consuming jobs have completed, unless temporary-file cleanup is disabled.

!!! note "Note"
    Check the process exit status and logs as well as output files. Empty DE tables, enrichment
    skip statuses, placeholder plots, and CellChat status objects need interpretation before they
    are treated as successful scientific results.

</div>
