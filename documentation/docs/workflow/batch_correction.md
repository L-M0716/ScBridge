# Batch Correction & Integration

ScBridge processes the merged normalized object using the selected integration method and compares
embeddings before and after correction.

The analysis aims to reduce technical batch-associated variation while retaining interpretable cell populations. The before-and-after embeddings provide diagnostics, not an automatic measure of correction quality.

## Input requirements

- **Expression object:** `normalization/seurat_merged_normalized.rds`, containing normalized expression data.
- **Batch metadata:** A cell-metadata column named by `integration.batch_key`, with batch assignments for the cells.
- **Sample and group metadata:** Sample identifiers and biological-group labels retained for downstream analyses.

## Processing workflow

1. **Initial assessment:** Inspect the available assay and variable-feature set, then calculate pre-correction PCA and UMAP.
2. **Batch correction:** Apply the selected integration method when multiple batches are present.
3. **Embedding assessment:** Calculate the corresponding post-correction embedding and compare batch distributions.
4. **Object export:** Save the processed Seurat object for variable-feature assessment and clustering.

## Integration methods

| Method value | Implementation | Representation |
|---|---|---|
| `seurat` | SCT preparation, anchor detection, and Seurat integration | Integrated assay and PCA |
| `harmony` | Harmony correction of PCA embeddings | Harmony reduction |
| `combat` | `sva::ComBat_seq` on RNA counts, followed by normalization and PCA | ComBat assay and PCA |

For datasets containing a single batch, no between-batch correction is applied; the pre-correction
PCA and UMAP embeddings are retained. Batch metadata remain required for diagnostic visualization,
including single-batch datasets.

## Parameter selection {#parameters}

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `integration.method` | `seurat` | Selection of the correction method |
| `integration.batch_key` | `batch` | Identification of batch assignments |
| `integration.n_pcs` | `30` | Dimensionality of PCA and subsequent embeddings |
| `integration.n_features` | `3000` | Variable-feature selection when no feature set is available |
| `normalization.seed` | `42` | Random-number initialization |

</div>

Seurat integration normally uses the existing variable-feature list. Changing `n_features`
does not necessarily replace that list; consensus-feature selection is described in [Normalization](normalization.md).

### Parameter guidance

- **Batch definition:** Set `integration.batch_key` to the technical grouping intended for correction. A condition label should not be substituted without considering the study design.
- **Method comparison:** If comparing integration methods or dimension counts, keep the input cells and upstream normalization consistent and use separate result directories. Compare the resulting structure as well as batch mixing.

## Run the analysis {#run-the-stage}

Integration runs as a dependency of `-t clustering` and later targets. There is no
`-t batch_correction` option in the launcher.

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
    -t clustering
```

## Outputs

Files are written under `RESULTS_DIR/batch_correction/`.

| File | Contents |
|---|---|
| `batch_corrected.rds` | Object passed to feature selection |
| <code>batch_<wbr>correction_<wbr>umap.<wbr>pdf</code> | Before/after UMAP by batch and cell-cycle phase when available |
| `batch_correction_pca.pdf` | Pre-correction PCA by batch |

## Inspect the results

<figure class="result-figure" markdown="1">

[![UMAP before and after correction, colored by batch and cell-cycle phase.](../assets/examples/integration.png)](../assets/examples/integration.png)

<figcaption markdown="1">

**Example report.** Compare batch representation within corresponding populations. The two embeddings have different coordinates, so distances and orientation are not directly comparable. Assess preservation of biological structure with marker evidence rather than batch mixing alone.


</figcaption>
</figure>

!!! note "Note"
    - **Inspect:** Compare the before-and-after batch-colored views in `batch_correction_umap.pdf`; inspect cell-cycle-phase panels when available. `batch_correction_pca.pdf` contains the pre-correction PCA view, not an elbow plot.
    - **Decide:** Examine whether differences between comparable populations persist and whether distinct populations collapse together. Use downstream marker evidence to assess biological preservation; the batch-colored view alone cannot establish this.
    - **Technical variation:** Compare pre- and post-correction embeddings using batch labels.
      Improved batch mixing should be assessed within comparable cell populations.
    - **Biological structure:** Evaluate the preservation of cell populations and marker-expression
      patterns alongside mixing. Greater overlap alone does not demonstrate improved integration.
    - **Experimental design:** Batch assignments represent technical groupings, whereas condition
      labels represent the biological comparison. When the two are confounded, correction may
      also attenuate variation relevant to the study.

Next: [Feature Selection](feature_selection.md).
