# Cell-Cell Communication

ScBridge uses CellChat to infer ligand-receptor communication potential between annotated cell groups.

This analysis asks which annotated cell groups show expression-based support for ligand-receptor communication. It evaluates communication potential rather than directly observed signaling.

## Input requirements

- An annotated Seurat object with expression data and cell-type metadata.
- A valid `cell_type_column`; the default is `singler_labels`.
- At least two cell types overall for communication inference.
- Sufficient cells in each group for the selected filtering threshold.

The ScBridge annotation stage stores final fine labels in `singler_labels`. For an external RDS,
ensure that the metadata column selected by `cell_type_column` contains the cell-group labels
intended for communication analysis.

## Processing workflow {#inference-workflow}

1. **Identity and database selection:** Select cell identities and the human or mouse CellChat database according to `organism`.
2. **Interaction screening:** Apply the requested database subset and identify overexpressed signaling genes and interactions.
3. **Probability estimation:** Estimate communication probabilities with `computeCommunProb(type = "triMean")`.
4. **Network aggregation:** Filter communications using `min_cells`, summarize pathways, and aggregate networks.
5. **Result visualization:** Generate network, pathway, and outgoing and incoming signaling-pattern reports.

When metadata contain two or more `group` values, the script runs inference by group and combines
successful objects. Groups with fewer than 50 cells are skipped. With no multi-group split,
inference uses the complete object.

## Parameter selection {#parameters}

Select the cell-identity column and database subset for communication inference.
The remaining settings control cell-group filtering, pathway visualization, and the number
of communication patterns requested for sending and receiving cells.

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `cell_type_column` | `singler_labels` | Cell-group identities |
| `db` | `all` | Full database or a CellChat-supported subset |
| `min_cells` | `10` | Minimum cell-group size for communication filtering |
| `pathways_of_interest` | `[]` | Optional pathway selection for visualization |
| `top_n_pathways` | `20` | Pathways displayed in summaries |
| `n_patterns_outgoing` | `4` | Requested outgoing pattern count |
| `n_patterns_incoming` | `4` | Requested incoming pattern count |

</div>

`[]` denotes an empty list: no specific pathways are requested for visualization. The report
then uses the first `top_n_pathways` entries in the inferred pathway list.

`min_cells` sets the minimum cell-group size used to filter inferred communications.
`pathways_of_interest` selects from the detected pathways; if no list is provided or none of its
entries match, the report uses the first `top_n_pathways` entries in the inferred pathway list.
`n_patterns_outgoing` and `n_patterns_incoming` specify the requested numbers of patterns
for sending and receiving cells, respectively.

### Parameter guidance

- **Cell identities:** Use `cell_type_column` to define the sender and receiver groups at the intended annotation resolution. Check cell-group sizes before interpreting effects of `min_cells`.
- **Reports:** `pathways_of_interest` and `top_n_pathways` affect pathway visualization; the requested incoming and outgoing pattern counts affect pattern analysis. These settings do not change the cell-type assignments.

## Run the analysis

This example resumes from an annotated object. Run from the repository root, adjust the paths,
and append `-n` to preview the planned jobs.

```bash
apptainer exec \
  -B "$PWD":/opt/scRNA_workflow \
  scRNA_seq.sif \
  bash /opt/scRNA_workflow/run_workflow \
    -I rds \
    -G annotation \
    -P /opt/scRNA_workflow/rds/seurat_final_annotated.rds \
    -C /opt/scRNA_workflow/config/config.yaml \
    -S /opt/scRNA_workflow/config/samples.tsv \
    -M /opt/scRNA_workflow/config/markers.tsv \
    -R /opt/scRNA_workflow/results_downstream \
    -t cellchat
```

## Outputs

Files are written under `RESULTS_DIR/cellchat/`.

| Files | Contents |
|---|---|
| `cellchat_object.rds` | Inference object, merged object, or skip/error status |
| `cellchat_final.rds` | Object from pattern analysis |
| `spreadsheets/cell_communication_network.csv` | Inferred communication network |
| `spreadsheets/LR_pairs_all.csv` | Extracted ligand-receptor pairs |
| `spreadsheets/pathway_summary.csv` | Pathway-level summary |
| `spreadsheets/centrality_scores.csv` | Network centrality scores |
| `plots/network/` | Interaction circle and heatmaps |
| `plots/pathways/` | Pathway and signaling-role plots |
| `plots/patterns/` | Outgoing and incoming signaling-pattern plots |

## Inspect the results

[![Cell-group interaction networks from the CellChat report.](../assets/examples/cellchat.png)](../assets/examples/cellchat.png)

**Example report.** Use the network overview to identify cell-group pairs for inspection in the interaction tables. The figure alone does not establish physical contact, causal signaling, or a statistically tested difference between conditions.


!!! note "Note"
    - **Inspect:** Compare the network views in `plots/network/` with `spreadsheets/cell_communication_network.csv`, then use `plots/pathways/` and `spreadsheets/LR_pairs_all.csv` to examine the contributing pathways and ligand-receptor pairs.
    - **Decide:** Check group composition and retained cell counts before attributing differences to signaling. Interpret `plots/patterns/` as summaries of shared sending and receiving patterns, and confirm that inference ran successfully for the groups being compared.
    - **Interpretation:** Inferred interactions represent expression-based communication
      potential, not direct measurements of signaling between cells.
    - **Assessment:** Evaluate cell-type assignments, group composition, and cell counts
      alongside network strength and the ligand-receptor pairs supporting each pathway.
    - **Limitations:** Skipped or unsuccessful analyses may produce status objects or message
      tables. Inspect the stage logs and output contents; file existence alone does not
      confirm successful inference, and missing interactions do not establish absent signaling.
