# Cell Annotation

ScBridge combines SingleR reference-based broad labels with lineage-restricted UCell marker scoring.
Final labels are assigned at the cluster level and mapped back to cells.

The broad-to-fine hierarchy makes the evidence for subtype assignments explicit. Its purpose is to refine supported identities within eligible lineages, not to force every cluster into the most specific available label.

## Input requirements

- A clustered Seurat object with `seurat_clusters`, RNA or SCT expression, and UMAP for reporting.
- `organism` set to `human` or `mouse`.
- A marker table supplied with `-M`; gene symbols must match the object.

## Marker table

| Column | Contents |
|---|---|
| `species` | Species tag, such as `Hs` or `Mm` |
| `main_cell_type` | Parent cell type |
| `sub_cell_type` | Subtype label |
| `marker_gene` | Marker gene symbol |

The reader also recognizes `official gene symbol` and `cell type` as aliases for
`marker_gene` and `main_cell_type`. Without subtype labels, the table cannot drive subtype refinement.
See [Input Data Overview](../inputs/index.md#marker-tables) for species handling.

## Processing workflow {#annotation-workflow}

1. **Reference-based labeling:** Apply SingleR to the clustered Seurat object to obtain per-cell reference labels and scores.
2. **Broad assignment:** Summarize per-cell labels and scores within each cluster to assign a broad cell type.
3. **Lineage matching:** Match broad labels to the parent cell types represented in the marker table.
4. **Subtype scoring:** Compare cluster-level median UCell scores among eligible subtypes within the matched lineage.
5. **Label assignment:** Assign supported subtype labels, retain broad labels when subtype evidence is insufficient, and map the final labels back to cells.

### Broad assignment

SingleR uses celldex `HumanPrimaryCellAtlasData` for human data or `MouseRNAseqData` for mouse data.
Pruned per-cell labels are summarized by cluster: the modal label and median positive score determine
the broad assignment. Clusters below the score threshold are labeled `Unknown`.

### Subtype refinement

Broad labels are matched to marker-table lineages using the script's mapping and matching rules.
Only subtypes belonging to the matched lineage compete. Median UCell scores are compared within
each cluster; insufficient evidence retains the broad label instead of forcing a subtype.
An `Unknown` broad assignment remains `Unknown`.

## Parameter selection

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `organism` | `human` | Reference species; supported values are `human` and `mouse` |
| <code>differential.<wbr>logfc_<wbr>threshold</code> | `0.25` | Log-fold-change threshold for marker testing in the annotation report |
| `differential.min_pct` | `0.1` | Minimum detected-cell fraction in either comparison group for marker testing |

</div>

These reporting parameters do not alter the SingleR or UCell assignment thresholds.

### Fixed assignment thresholds {#evidence-settings}

<div class="parameter-table" markdown="1">

| Parameter | Fixed value | Description |
|---|---|---|
| `HPCA_CONF_THRESH` | `0.20` | Cluster-level broad assignment |
| `UCELL_MINMARKERS` | `2` | UCell signature eligibility |
| `UCELL_MAXRANK` | `1500` | Signature scoring |
| `UCELL_MARGIN` | `0.03` | Subtype acceptance, with a positive best score |

</div>

These are fixed constants in the annotation script, not user-configurable parameters.

### Parameter guidance

- **Reference and markers:** Match `organism`, gene identifiers, and marker-table species before interpreting an assignment. Use subtype markers appropriate to the sampled tissue and the candidate parent lineage.
- **Reporting thresholds:** Adjust `differential.logfc_threshold` and `differential.min_pct` only to change marker reporting. These values do not relax the fixed SingleR or UCell acceptance rules.

## Run the analysis {#run-the-stage}

- **Run location:** Repository root, with `scRNA_seq.sif` available there.
- **Paths:** Replace the example paths with container-accessible paths for your files.
- **Preview:** Append `-n` to preview jobs; setup may still generate runtime files.

```bash
apptainer exec \
  -B "$PWD":/opt/scRNA_workflow \
  scRNA_seq.sif \
  bash /opt/scRNA_workflow/run_workflow \
    -I matrix \
    -D /opt/scRNA_workflow/matrix_files \
    -C /opt/scRNA_workflow/config/config.yaml \
    -S /opt/scRNA_workflow/config/samples.tsv \
    -M /opt/scRNA_workflow/config/markers.tsv \
    -R /opt/scRNA_workflow/results_matrix \
    -t annotation
```

## Outputs

Files are written under `RESULTS_DIR/annotation/`.

| File | Contents |
|---|---|
| <code>seurat_<wbr>final_<wbr>annotated.<wbr>rds</code> | Final object with `celltype_main` and `celltype_fine` |
| <code>spreadsheets/<wbr>auto_<wbr>annotation_<wbr>scores.<wbr>csv</code> | Cluster annotation evidence and scores |
| <code>reports/<wbr>auto_<wbr>annotation_<wbr>report.<wbr>pdf</code> | Broad/fine annotation diagnostics |
| <code>reports/<wbr>umap_<wbr>final_<wbr>annotated.<wbr>pdf</code> | Annotated UMAP |
| <code>reports/<wbr>celltype_<wbr>composition.<wbr>pdf</code> | Cell-type composition plots |
| <code>spreadsheets/<wbr>celltype_<wbr>composition.<wbr>csv</code> | Composition table |
| <code>reports/<wbr>marker_<wbr>gene_<wbr>heatmap.<wbr>pdf</code> | Marker-expression heatmap |

The annotation script also sets `singler_labels` to the fine labels for downstream compatibility.

## Inspect the results

[![Original clusters, major cell types, and cell subtypes on the same UMAP.](../assets/examples/annotation.png)](../assets/examples/annotation.png)

**Example report.** Compare the same regions across the three panels. Broad monocyte labels are refined into classical and non-classical monocytes, while some labels remain broad or Unknown. Review marker expression and cluster scores before accepting these assignments.


!!! note "Note"
    - **Broad labels:** In `reports/auto_annotation_report.pdf`, compare the broad-label view with the cluster-level evidence in `spreadsheets/auto_annotation_scores.csv`.
    - **Subtype refinement:** Compare broad and fine labels and identify which assignments became more specific. A retained broad label or `Unknown` is an allowed outcome when support is insufficient.
    - **Expression evidence:** Review `reports/marker_gene_heatmap.pdf` alongside the score table and marker definitions; a single expressed marker is not sufficient evidence for a subtype.
    - **Sample context:** Use `reports/umap_final_annotated.pdf` and the composition outputs to assess label distribution before selecting populations for downstream analysis.
    - **Interpretation:** Labels are computational assignments, not definitive biological identities. Review score support, multiple markers, cluster heterogeneity, sample composition, and tissue context. Reference coverage and marker quality limit subtype resolution; do not relabel an uncertain cluster solely to eliminate `Unknown`.

Continue with [Differential Expression](../downstream/differential.md), [Trajectory Inference](../downstream/monocle3.md), or [Cell-Cell Communication](../downstream/cellchat.md). [Functional Enrichment](../downstream/enrichment.md) uses the differential-expression results.
