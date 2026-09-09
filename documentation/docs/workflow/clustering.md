# Clustering

ScBridge constructs a neighbor graph, assigns clusters, and generates a UMAP embedding from the
processed Seurat object.

Clusters define the groups summarized during cell annotation. A useful clustering result should support coherent expression patterns while allowing sample and batch composition to be examined.

## Input requirements

- `feature_selection/seurat_feature.rds`, with suitable expression data or reductions.
- Sample metadata containing `sample_id` and `group`.
- Matching sample identifiers in the object and metadata; a `batch` column for preview plots.

## Processing workflow {#processing-steps}

1. **Metadata assignment:** Map biological groups from the sample table to cells.
2. **Reduction selection:** Use Harmony dimensions when the selected method is Harmony and its reduction exists; otherwise use PCA.
3. **PCA preparation:** In the PCA branch, calculate principal components if the reduction is missing or contains fewer dimensions than requested.
4. **Graph construction and clustering:** Build neighbors using `FindNeighbors`, then assign clusters using `FindClusters`.
5. **Embedding visualization:** Generate UMAP from the same reduction and export cluster, batch, and group views.

## Parameter selection {#parameters}

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `integration.n_pcs` | `30` | Dimensions used by neighbors and UMAP |
| `clustering.n_neighbors` | `20` | Neighbor-graph size |
| <code>clustering.<wbr>resolution.<wbr>default</code> | `0.8` | Clustering resolution |
| `clustering.algorithm` | `Leiden` | `Leiden` or `Louvain` |
| `clustering.umap.min_dist` | `0.3` | UMAP local compactness |
| `clustering.umap.spread` | `1.0` | UMAP embedding spread |
| `normalization.seed` | `42` | Random seed |

</div>

Use the documented algorithm names: an unrecognized value currently falls back to Louvain.

### Parameter guidance

- **Graph structure:** `integration.n_pcs` and `clustering.n_neighbors` affect the neighbor graph. Assess their effects separately from changes to clustering resolution.
- **Resolution:** Compare plausible resolutions using consistent input cells, reduction, and seed settings in separate result directories. Review markers and sample representation rather than choosing the most visually separated UMAP.
- **UMAP display:** `clustering.umap.min_dist` and `clustering.umap.spread` affect the embedding, not the preceding `FindClusters` assignments.

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
    -t clustering
```

## Outputs

Files are written under `RESULTS_DIR/clustered/`.

| File | Contents |
|---|---|
| `seurat_clustered.rds` | Seurat object with `seurat_clusters`, group metadata, and UMAP |
| <code>clustering_<wbr>preview_<wbr>umap.<wbr>pdf</code> | Cluster/batch overview and per-group cluster plots |

## Inspect the results

<figure class="result-figure" markdown="1">

[![The same UMAP colored by cluster and batch.](../assets/examples/clustering.png)](../assets/examples/clustering.png)

<figcaption markdown="1">

**Example report.** Compare cluster boundaries with batch representation, then inspect marker expression before assigning biological labels. Separated UMAP regions and numbered clusters are not, by themselves, evidence of distinct cell types.


</figcaption>
</figure>

!!! note "Note"
    - **Inspect:** Open `clustering_preview_umap.pdf` for the cluster and batch views and the per-group cluster panels. Identify clusters concentrated in one batch or group and check their cell counts and expression evidence.
    - **Decide:** Use marker evidence from the annotation stage to assess whether apparent subdivisions are interpretable. If revising the graph or resolution, repeat the affected downstream stages so that annotation matches the new clusters.
    - **Interpretation:** Compare resolutions using marker coherence, sample representation, and stability of major populations. Higher resolution is not automatically more informative. UMAP separation alone does not establish a distinct cell type, and cluster numbers have no fixed biological meaning.

Next: [Cell Annotation](annotation.md).
