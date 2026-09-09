# Normalization

ScBridge performs sample-wise normalization using SCTransform v2, followed by cross-sample
variable-feature selection and merging of normalized Seurat objects for downstream integration.

The resulting SCT assays and consensus feature set provide a common starting point for cross-sample integration. Normalization and batch correction address different parts of the workflow.

## Input requirements

- **Expression data:** Filtered Seurat objects containing an RNA counts layer and cell-level metadata.
- **Cell-cycle reference:** A species-matched RDS containing `s.genes` and `g2m.genes` gene sets.
- **Covariates:** Cell-metadata columns corresponding to the variables in `normalization.vars_to_regress`.

## Processing workflow {#normalization-procedure}

1. **RNA log-normalization:** Generate a log-normalized RNA data layer using a scale factor of 10,000.
2. **Cell-cycle scoring:** Calculate S-phase and G2M-phase scores when at least five genes from each reference set are present in the expression data.
3. **SCTransform:** Apply SCTransform v2 with `glmGamPoi` and the specified regression covariates to each sample.
4. **Cross-sample feature selection:** Use `SelectIntegrationFeatures` to select consensus variable features, then merge the normalized objects and assign the selected feature set.

RNA log-normalization provides the expression values used for cell-cycle scoring, whereas
SCTransform generates the SCT assay used in subsequent processing. The two transformations are
performed within the same normalization procedure.

## Parameter selection {#parameters}

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| <code>normalization.<wbr>vars_<wbr>to_<wbr>regress</code> | `["percent.mt"]` | Metadata variables regressed during SCTransform |
| <code>normalization.<wbr>n_<wbr>variable_<wbr>genes</code> | `3000` | Consensus features selected across normalized samples |
| `normalization.seed` | `42` | Random seed |
| <code>qc.<wbr>gene_<wbr>filtering.<wbr>min_<wbr>cells</code> | `7` | Minimum cell detection count passed to SCTransform |
| `input.cell_cycle_genes` | Species-matched RDS | S and G2M gene signatures |

</div>

!!! note "Feature selection and covariate adjustment"
    - **Feature sets:** Per-sample SCTransform specifies `variable.features.n = 3000` and
      `return.only.var.genes = FALSE`; the SCT output is therefore not restricted to the selected
      variable features. The cross-sample consensus-feature count is controlled separately by
      `normalization.n_variable_genes`.
    - **Covariate adjustment:** Regression is restricted to the metadata variables listed in
      `normalization.vars_to_regress`. Calculating cell-cycle scores does not, by itself, remove
      cell-cycle-associated variation.

### Parameter guidance

- **Regression covariates:** Keep `normalization.vars_to_regress` tied to an explicit technical concern and ensure the named metadata exist. Adding S and G2M scores changes the regression model; calculating those scores alone does not regress them.
- **Consensus features:** Change `normalization.n_variable_genes` when evaluating the sensitivity of downstream structure to feature selection. This does not change the script's per-sample 3,000-feature setting.

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
    -t normalization
```

## Outputs

Files are written under `RESULTS_DIR/normalization/`.

| File | Contents |
|---|---|
| <code>seurat_<wbr>merged_<wbr>normalized.<wbr>rds</code> | Merged Seurat object with RNA and SCT assays |
| `consensus_hvg_list.txt` | Consensus variable-feature list |
| <code>normalization_<wbr>qc_<wbr>figures.<wbr>pdf</code> | Normalized count distribution by sample |
| <code>&lt;sample_<wbr>id&gt;/<wbr>seurat_<wbr>normalized.<wbr>rds</code> | Temporary per-sample intermediate; may be removed by Snakemake |

## Inspect the results

<figure class="result-figure" markdown="1">

[![Normalized count distributions across five samples.](../assets/examples/normalization.png)](../assets/examples/normalization.png)

<figcaption markdown="1">

**Example report.** Compare sample distributions and investigate unexpected differences alongside the assay and normalization logs. Similar distributions are not required, and this plot alone does not demonstrate removal of batch effects.


</figcaption>
</figure>

!!! note "Note"
    - **Inspect:** Review `normalization_qc_figures.pdf` for sample-level count distributions and `consensus_hvg_list.txt` for the selected feature set. Confirm that the merged object contains the RNA and SCT assays expected by integration.
    - **Decide:** Resolve absent covariates, unexpected feature loss, or unsuccessful normalization before integration. Distribution differences alone do not justify removing a biological covariate.
    - **Sample comparability:** Assess normalized count distributions and the representation of
      selected features across samples. Normalization alone does not establish that technical
      batch effects have been removed.
    - **Biological variation:** Covariate adjustment should be justified by the experimental design.
      Variables associated with both technical effects and the biological process of interest may
      remove relevant signal when regressed.
    - **Diagnostic availability:** Interpretation should be based on the generated plots and
      normalization logs. When plot generation is unsuccessful, the report may contain a status
      message rather than a diagnostic figure.

Next: [Batch Correction & Integration](batch_correction.md).
