# Differential Expression

ScBridge compares biological groups within each annotated cell type using sample-level pseudobulk
counts and DESeq2.

This analysis addresses which genes differ between the selected biological groups within each cell type. The unit of replication is the sample, not the individual cell.

## Analysis design

Counts are summed across cells belonging to the same sample and cell type. DESeq2 uses the
grouping column as the design term and compares one non-control group against the configured control.

- **Replicates:** At least two samples per group for each tested cell type.
- **Cell coverage:** At least `min_cells_group` cells per group within that cell type.
- **Counts:** RNA counts are preferred; SCT counts are used only when RNA is absent.
- **Labels:** `celltype_fine`, then `singler_labels`, then `manual_annotation`.
- **Groups:** Metadata must map sample identifiers to the requested comparison groups.

!!! note "Note"
    The current script runs one contrast. With more than two groups, it selects the first
    non-control group in the sample table. Use an explicit two-group study subset when preparing
    a comparison. The design is group-only: paired subjects, batch covariates, and interaction
    terms are not modeled by this wrapper.

## Parameter selection {#testing-and-thresholds}

Genes must have at least one count in at least the smaller group's sample count.
DESeq2 receives `alpha = 0.05` and `lfcThreshold = differential.logfc_threshold`.
Results with missing adjusted P values are removed; `significant` marks adjusted P values below 0.05.

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| <code>differential.<wbr>group_<wbr>by_<wbr>column</code> | `group` | Group metadata column |
| <code>differential.<wbr>control_<wbr>group</code> | `Adult peripheral blood` | Reference-group label matching the sample metadata |
| <code>differential.<wbr>logfc_<wbr>threshold</code> | `0.25` | DESeq2 log2 fold-change threshold |
| <code>differential.<wbr>min_<wbr>cells_<wbr>group</code> | `20` | Minimum cells per group and cell type |
| <code>differential.<wbr>visualization.<wbr>adj_<wbr>pval_<wbr>cutoff</code> | `0.05` | Plot significance cutoff |
| <code>differential.<wbr>visualization.<wbr>logfc_<wbr>threshold</code> | `0.5` | Plot effect-size cutoff |
| <code>differential.<wbr>visualization.<wbr>top_<wbr>n</code> | `20` | Genes selected for summaries |

</div>

`differential.min_pct` is used by annotation marker reporting, not by this DESeq2 test.
There is no configurable `de_method: wilcox` route in the current differential script.

### Parameter guidance

- **Contrast definition:** Confirm that `differential.control_group` exactly matches the intended reference label and that each tested cell type has the required sample replication.
- **Testing versus display:** Distinguish the DESeq2 effect-size threshold from the visualization cutoffs. Changing a volcano-plot cutoff changes the highlighted results, not the fitted model.

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
    -t differential
```

## Outputs

Files are written under `RESULTS_DIR/differential/`.

| Files | Contents |
|---|---|
| <code>spreadsheets/<wbr>differential_<wbr>expression_<wbr>results.<wbr>csv</code> | Gene-level results by cell type and comparison |
| <code>spreadsheets/<wbr>top_<wbr>de_<wbr>genes_<wbr>summary.<wbr>csv</code> | Selected top genes |
| <code>reports/<wbr>de_<wbr>volcano_<wbr>plots.<wbr>pdf</code> | Effect-size and significance plots |
| `reports/de_heatmap.pdf` | Differential-expression heatmap |

`avg_log2FC` stores the DESeq2 log2 fold change, positive for the non-control group.
`p_val_adj` stores the adjusted P value. The primary table is not restricted to significant rows.

## Inspect the results

[![CD4 T-cell volcano plot for umbilical cord blood versus adult peripheral blood.](../assets/examples/differential.png)](../assets/examples/differential.png)

**Example report.** The horizontal axis shows log2 fold change and the vertical axis shows the negative log10 adjusted P value. The report highlights genes passing its displayed thresholds. Check the underlying counts and replication, especially for extreme fold changes.


!!! note "Note"
    - **Inspect:** Start with `spreadsheets/differential_expression_results.csv`, checking the comparison and cell-type columns, effect direction, and adjusted P values. Use `reports/de_volcano_plots.pdf` and `reports/de_heatmap.pdf` to examine the corresponding expression patterns.
    - **Decide:** Confirm replication and contrast direction before selecting genes for enrichment. A heatmap or a large fold change does not replace the sample-level statistical evidence.
    - **Interpretation:** Results quantify group-associated expression differences within each
      tested cell type. Interpret effect sizes together with adjusted P values.
    - **Assessment:** Evaluate sample-level replication and cell coverage for each comparison.
      Sample identifiers must represent independent biological replicates; additional cells
      do not substitute for additional replicates.
    - **Limitations:** Empty results may reflect insufficient samples, insufficient cells, or
      unsuccessful testing. Inspect the per-cell-type log before interpreting an empty table
      as evidence of no differential expression.

Next: [Functional Enrichment](enrichment.md) for GO and KEGG over-representation analysis.
