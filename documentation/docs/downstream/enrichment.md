# Functional Enrichment

ScBridge performs GO and offline KEGG over-representation analysis on differential-expression
results, separately for each comparison and cell type.

The biological question is whether selected differential-expression genes are over-represented in annotated processes or pathways. Upregulated and downregulated sets are evaluated separately as well as together.

## Gene selection

1. Read `differential/spreadsheets/differential_expression_results.csv`.
2. Select genes with `p_val_adj < pval_cutoff` and `abs(avg_log2FC) > logfc_threshold`.
3. Form three gene sets: all selected genes, upregulated genes, and downregulated genes.
4. Convert symbols to unique Entrez IDs using the organism-specific annotation database.
5. Test eligible sets and export tables and figures.

Sets with fewer than five genes, or fewer than five mapped Entrez IDs, are skipped and recorded
in the analysis summary.

## GO and KEGG

| Analysis | Implementation | Resources |
|---|---|---|
| GO biological process | `enrichGO` | Organism-specific OrgDb |
| GO molecular function | `enrichGO` | Organism-specific OrgDb |
| GO cellular component | `enrichGO` | Organism-specific OrgDb |
| KEGG pathways | `enricher` | Offline `term2gene` and `term2name` mappings |

Human data use `org.Hs.eg.db` and human KEGG mappings; mouse data use `org.Mm.eg.db` and mouse mappings.
Both enrichment calls use Benjamini-Hochberg adjustment.

!!! note "Note"
    The current calls do not supply a custom `universe`. The background therefore follows the
    annotation mappings used by the enrichment functions, not a user-defined list of detected
    or tested genes. Consider this when interpreting results from strongly filtered datasets.

## Parameter selection {#parameters}

Gene-selection thresholds determine which differential-expression results enter enrichment.
Term-size limits constrain the tested gene sets, while the enrichment cutoffs control result
retention. `top_n_terms` affects visualization only; it does not change the enrichment tests.

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `pval_cutoff` | `0.05` | Input DEG adjusted-P filter and enrichment P-value cutoff |
| `qval_cutoff` | `0.1` | Enrichment q-value cutoff |
| `logfc_threshold` | `0.5` | Absolute input log2 fold-change threshold |
| `min_gene_size` | `10` | Minimum term size |
| `max_gene_size` | `500` | Maximum term size |
| `top_n_terms` | `20` | Number of enrichment terms selected for visualization |

</div>

### Parameter guidance

- **Gene selection:** Assess `pval_cutoff` and `logfc_threshold` together with the number of mapped genes. An extremely small set can be skipped rather than yield an informative negative result.
- **Term selection:** Interpret term-size and q-value filters as part of the analysis definition. Changing `top_n_terms` changes the number displayed, not the underlying enrichment result.

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
    -t enrichment
```

Enrichment depends on the differential-expression result table. If that table is missing,
Snakemake schedules the upstream differential-expression step before enrichment.

## Outputs

Files are written under `RESULTS_DIR/enrichment/`.

| Files | Contents |
|---|---|
| `GO/` and `KEGG/` | Tables by comparison, cell type, gene set, and ontology |
| `enrichment_summary.csv` | Gene counts, term counts, and skip statuses |
| `enrichment_summary.txt` | Human-readable summary |
| `enrichment_results_all.rds` | Collected enrichment objects |
| `unmapped_genes.csv` | Unmapped gene symbols |
| `GO_visualization/` and `KEGG_visualization/` | Enrichment plots |

## Inspect the results

[![GO biological-process enrichment for the classical-monocyte all-gene selection.](../assets/examples/enrichment.png)](../assets/examples/enrichment.png)

**Example report.** Read GeneRatio together with the contributing gene count and adjusted P value. This all-gene selection combines expression-change directions; enriched terms do not establish pathway activation or inhibition.


!!! note "Note"
    - **Inspect:** Read `enrichment_summary.csv` before opening the GO and KEGG plots. Check gene-set size, mapped identifiers, and whether each comparison was tested or skipped; then inspect the genes contributing to the reported terms.
    - **Decide:** Interpret related terms together and separate upregulated from downregulated sets. Revisit identifier mapping or resource availability before changing cutoffs to obtain more terms.
    - **Interpretation:** Results identify functional terms over-represented among the selected
      genes relative to the annotation background. Enrichment does not establish pathway activation.
    - **Assessment:** Evaluate gene-set size, identifier-mapping success, selection thresholds,
      and the genes contributing to each term. Confirm that the species-specific resources
      are appropriate for the input.
    - **Limitations:** Missing offline KEGG resources or insufficient mapped genes can cause
      analyses to be skipped. Inspect the analysis summary and logs; an empty table does not
      establish the absence of functional enrichment.
