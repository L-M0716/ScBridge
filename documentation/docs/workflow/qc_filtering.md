# Quality Control & Filtering

The QC path calculates per-cell quality metrics, detects doublets, and removes cells that fail configured
thresholds before normalization.

This stage establishes which cells and genes enter the shared analysis. Retention decisions should be assessed per sample before interpreting differences between biological groups.

## Input requirements

- Gene-by-cell count matrices with matching barcode and feature identifiers.
- Sample identifiers supplied through `-S`, matching the input sample directories.
- Mitochondrial and gene-exclusion patterns appropriate for the gene identifiers in the data.

## Processing workflow {#analysis-workflow}

1. **Quality assessment:** Calculate detected gene counts, total counts, and mitochondrial percentages from each sample's count matrix.
2. **Doublet detection:** Run scDblFinder separately for each sample and add doublet predictions to cell metadata.
3. **Cell filtering:** Retain cells that pass the QC thresholds and are not classified as doublets.
4. **Gene filtering:** Apply the minimum cell-detection threshold and remove genes matching the exclusion patterns.
5. **Object export:** Save the filtered Seurat object and filtering summaries for each sample.

### QC metrics

| Metric | Meaning |
|---|---|
| `nFeature_RNA` | Detected genes per cell |
| `nCount_RNA` | Total counts per cell |
| `percent.mt` | Percentage assigned to genes matching the mitochondrial pattern |
| `qc_status` | Combined pass/fail result for the four thresholds |

### Doublet detection

ScBridge runs scDblFinder separately for each sample. The expected doublet rate and random seed
are controlled by `qc.expected_doublet_rate` and `normalization.seed`, with default values of
`0.06` (6%) and `42`, respectively. The expected rate is a model input, not a fixed proportion of
cells that will be removed.

Results are stored in `doublet_score`, `doublet_class`, and `doublet_call` within cell metadata;
`doublet_call` is `1` for predicted doublets and `0` otherwise.

`qc.doublet_method` records the method label; changing it does not switch the implemented detector.

## Parameter selection {#qc-parameters}

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `qc.min_genes` | `200` | Minimum detected genes per cell |
| `qc.max_genes` | `5000` | Maximum detected genes per cell |
| `qc.min_counts` | `700` | Minimum total counts per cell |
| `qc.max_mito` | `10` | Maximum mitochondrial count percentage per cell (%) |
| `qc.mito_pattern` | `^MT-` | Pattern used to identify mitochondrial genes |
| <code>qc.<wbr>expected_<wbr>doublet_<wbr>rate</code> | `0.06` | Expected doublet rate supplied to scDblFinder |
| <code>qc.<wbr>gene_<wbr>filtering.<wbr>min_<wbr>cells</code> | `7` | Minimum number of cells in which a retained gene is detected |
| <code>qc.<wbr>gene_<wbr>filtering.<wbr>exclude_<wbr>genes</code> | `^Rp[sl]`, `^Mt-`, `Malat1` | Gene-exclusion patterns applied during filtering |
| `normalization.seed` | `42` | Random seed for doublet detection |

</div>

Cells pass QC when the four cell-level thresholds (`min_genes`, `max_genes`, `min_counts`, and
`max_mito`) are satisfied, including the boundary values. Filtering retains
cells with `qc_status = pass` and `doublet_call = 0`.

Gene filtering then retains genes detected in at least `qc.gene_filtering.min_cells` cells
(`7` in the supplied configuration) and removes genes matching `qc.gene_filtering.exclude_genes`
(`^Rp[sl]`, `^Mt-`, and `Malat1`).

!!! important "QC considerations"
    - **Dataset-specific thresholds:** Review per-sample gene counts, total counts, and mitochondrial
      percentages before choosing cutoffs. The configured values are starting settings, not universal thresholds.
    - **Gene-name matching:** Confirm that `qc.mito_pattern` and gene-exclusion patterns match the
      identifiers and capitalization in the input. The supplied mitochondrial pattern is `^MT-`;
      mitochondrial matching is case-sensitive, whereas gene-exclusion matching is case-insensitive.
    - **Doublet calls:** Interpret predicted doublets alongside QC metrics and, when available,
      cluster structure and marker expression. A doublet prediction is not experimental confirmation.
    - **Retention checks:** Review the number of retained cells and genes for each sample before
      continuing to normalization, especially when filtering removes markedly different proportions across samples.

### Parameter guidance

- **Cell thresholds:** Use the per-sample gene-count, count-depth, and mitochondrial distributions to assess `qc.min_genes`, `qc.max_genes`, `qc.min_counts`, and `qc.max_mito`; do not transfer cutoffs between datasets without checking those distributions.
- **Gene identifiers:** Confirm the mitochondrial pattern before interpreting mitochondrial percentages. Changing a threshold does not correct an identifier-pattern mismatch.

## Run the analysis {#run-the-stage}

This example performs cell-level QC, doublet detection, and cell/gene filtering from count matrices,
stopping before normalization.

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
    -t filtering
```

The `-t filtering` target includes QC assessment, doublet detection, and cell and gene filtering.

## Outputs {#filtering-outputs}

Results are saved under `RESULTS_DIR/filtering/<sample_id>/`.

| File | Contents |
|---|---|
| `seurat_filtered.rds` | Filtered Seurat object for normalization or filtering-stage RDS input |
| `filtering_stats.csv` | Filtering summary statistics |
| `filtering_report.pdf` | Filtering diagnostic plots |

## Inspect the results

=== "Sample QC"

    <figure class="result-figure" markdown="1">

    [![Gene counts and mitochondrial percentages across five samples.](../assets/examples/qc.png)](../assets/examples/qc.png)

    <figcaption markdown="1">

    **Example report.** Compare the distributions sample by sample before choosing filtering thresholds. A shifted distribution can reflect cell composition or technical quality and needs further assessment.


    </figcaption>
    </figure>

=== "Doublets"

    <figure class="result-figure" markdown="1">

    [![Predicted singlets and doublets for GSM9101188.](../assets/examples/doublet.png)](../assets/examples/doublet.png)

    <figcaption markdown="1">

    **Example report.** Red points are predicted doublets; the report gives a 9.66% predicted doublet rate for GSM9101188. This is a model-derived classification, not an experimentally measured rate.


    </figcaption>
    </figure>

=== "Filtering"

    <figure class="result-figure" markdown="1">

    [![Retained and filtered cells in GSM9101188.](../assets/examples/filtering.png)](../assets/examples/filtering.png)

    <figcaption markdown="1">

    **Example report.** Check the retained and filtered distributions against the displayed cutoffs and the filtering statistics. This sample illustrates the report structure; its thresholds are not a universal recommendation.


    </figcaption>
    </figure>


!!! note "Note"
    - **Inspect:** Open each sample's `filtering_report.pdf` to compare retained and removed cells across gene-count and mitochondrial distributions and the count-versus-gene plot. Use `filtering_stats.csv` to quantify cell and gene retention.
    - **Decide:** Investigate disproportionate loss in individual samples before continuing. Review the relevant threshold, doublet calls, and gene-exclusion settings together; do not change thresholds solely to equalize retained cell counts.

Next: [Normalization](normalization.md).
