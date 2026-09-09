# Feature Selection

ScBridge examines variable features in the batch-corrected Seurat object and summarizes their
expression patterns before clustering. The stage uses available assay-specific feature sets,
updates variability statistics, and exports a feature report together with the processed object.

This stage summarizes the feature set carried forward from earlier processing. It supports inspection before clustering; it is not a separate interface for selecting an optimal number of principal components.

## Input requirements

- The batch-corrected object at `RESULTS_DIR/batch_correction/batch_corrected.rds`.
- An SCT or RNA assay for expression summaries; otherwise the original default assay is used.
- Existing UMAP coordinates for feature-expression maps. This stage does not calculate a new embedding.

## Processing workflow {#variable-feature-assessment}

1. **Assay selection:** Use SCT for expression summaries when available, otherwise RNA; if neither is present, retain the original default assay.
2. **Feature-set selection:** Use the integrated assay's variable features when available. Otherwise, use the visualization assay's existing feature list or estimate a list with `FindVariableFeatures` when none is available.
3. **Variability assessment:** Recompute variability statistics using the `vst` method with up to 3,000 features. The feature list selected in the preceding step remains the basis for the expression summaries.
4. **Expression visualization:** Summarize selected feature expression and, where suitable embeddings are available, display feature-expression maps.
5. **Object export:** Restore the original default assay and save the Seurat object for subsequent clustering.

## Implementation settings {#method-settings}

These settings are specified in the feature-selection script.

<div class="parameter-table" markdown="1">

| Setting | Value | Description |
|---|---|---|
| Initial feature-count fallback | `2000` | Used only when no variable-feature list is available |
| Variability-statistics feature count | `3000` | Upper limit, capped at the number of available genes |
| Variability estimation method | `vst` | Method used when updating variability statistics |

</div>

These settings are separate from the consensus-feature count used during normalization.

### Parameter guidance

- **Feature-count changes:** Adjust the cross-sample consensus count in [Normalization](normalization.md#parameters) when testing alternative feature sets. The fallback and visualization counts below are script settings, not additional command-line options.
- **Interpretation:** Inspect the represented genes and their expression patterns rather than assuming a larger feature set provides more biological information.

## Run the analysis {#run-the-stage}

The stage runs between integration and clustering when targeting `clustering` or a later analysis.
There is no `-t feature_selection` option.

See the container command in [Clustering](clustering.md#run-the-stage).

## Outputs

Files are written under `RESULTS_DIR/feature_selection/`.

| File | Contents |
|---|---|
| `seurat_feature.rds` | Seurat object passed to clustering |
| <code>hvg_<wbr>visualization_<wbr>report.<wbr>pdf</code> | Selected-feature expression summaries and available feature-expression maps |

## Inspect the results

<figure class="result-figure" markdown="1">

[![Selected feature values on the existing UMAP.](../assets/examples/features.png)](../assets/examples/features.png)

<figcaption markdown="1">

**Example report.** Review where each feature varies across the embedding. Color scales differ between panels and include negative values for some features; these are not raw molecule counts, and color intensity should not be compared directly across genes.


</figcaption>
</figure>

!!! note "Note"
    - **Inspect:** Open `hvg_visualization_report.pdf` for selected-feature expression summaries and maps on an existing embedding. Check the stage log if the report contains fewer panels than expected.
    - **Decide:** Investigate unexpected dominant genes in the input feature set and upstream filtering or normalization. This report does not provide a validated optimal HVG count or a PCA-dimension recommendation.
    - **Biological relevance:** Expression variability may reflect cell identity, cell state, or
      technical effects. Inclusion in a variable-feature set does not establish cell-type specificity.
    - **Expression context:** Interpret feature summaries alongside expression maps and sample
      composition. Marker-based cell annotation is performed in the subsequent annotation stage.
    - **Report scope:** Feature maps use existing embeddings rather than newly calculated PCA or
      UMAP coordinates. Available panels depend on the assays, features, and reductions in the
      object; consult the stage log when expected plots are absent.

Next: [Clustering](clustering.md).
