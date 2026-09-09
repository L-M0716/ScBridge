# Trajectory Inference

ScBridge uses Monocle3 to learn transcriptional trajectories and order annotated cells in pseudotime.

Trajectory inference asks how the selected cells vary along a connected transcriptional process. Population selection and the root reference define the scope and direction of the interpretation.

## Input and cell selection

The input is an annotated Seurat object with expression counts and `seurat_clusters`.
Cell labels preferentially use `celltype_fine`, with `singler_labels` as a fallback.
`clusters_to_use` accepts `all` or a comma-separated string of Seurat cluster IDs.

Choose populations relevant to the same biological process. Including every annotated cell type
is not necessarily appropriate for a single trajectory question.

## Processing workflow {#processing-steps}

1. **Cell selection:** Subset the requested Seurat clusters and construct a Monocle3 `cell_data_set` from expression counts and metadata.
2. **Preprocessing:** Apply size-factor normalization and dimensionality reduction with `preprocess_cds`, then calculate UMAP coordinates.
3. **Graph construction:** Cluster cells with `cluster_cells` and learn the principal graph with `learn_graph`.
4. **Pseudotime ordering:** Select root cells or a root principal node, then order cells with `order_cells`.
5. **Gene association:** Run `graph_test` on the principal graph and retain genes below the specified q-value threshold.
6. **Result export:** Save the trajectory object, analysis summary, selected gene statistics, and visualizations.

## Trajectory root specification {#root-selection}

The root defines the reference point for pseudotime ordering. A biologically justified cell-type
label or Seurat cluster can be supplied; otherwise, the script selects a root computationally.

<div class="parameter-table root-selection-table" markdown="1">

| Selection mode | Condition | Procedure |
|---|---|---|
| Cell-type reference | `root_cell_type` is specified and cell-type metadata are available | Use cells matching the requested annotation label |
| Cluster reference | No cell-type reference is applied and `root_cluster` is specified | Use cells matching the requested Seurat cluster ID |
| Automatic selection | Neither reference is specified, or the selected reference matches no cells | Select a principal-graph node using cell assignments; use the UMAP-based fallback if no valid node is available |

</div>

!!! note "Note"
    The cell-type reference takes precedence when its metadata are available. If the requested label
    matches no cells, the script proceeds to automatic selection rather than trying `root_cluster`.
    Automatic selection uses the principal node with the most assigned cells. If no valid node is
    available, it uses up to 500 cells from the cluster with the lowest median UMAP1 coordinate.
    These rules define a computational reference, not an experimentally established starting state.

## Parameter selection {#parameters}

Choose the cell population and root reference according to the biological question.
The preprocessing and neighborhood settings control trajectory construction, whereas
`q_value_threshold` filters graph-associated genes after the graph has been learned.

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `clusters_to_use` | `all` | Seurat clusters included in the analysis; accepts `all` or a comma-separated string of cluster IDs |
| `root_cell_type` | `""` | Optional annotation label identifying root cells |
| `root_cluster` | `""` | Optional Seurat cluster ID identifying root cells when no cell-type reference is applied |
| `num_dim` | `50` | Number of dimensions used by Monocle3 preprocessing |
| `resolution` | `0.001` | Resolution used by Monocle3 clustering before principal-graph learning |
| `k` | `20` | Number of nearest neighbors used by Monocle3 clustering |
| `q_value_threshold` | `0.05` | Maximum q value for retaining graph-associated genes; genes must fall strictly below this threshold |

</div>

`""` denotes an empty string: no root cell type or cluster is specified. When both root settings
are empty, the workflow selects a root automatically as described above.

### Parameter guidance

- **Population and root:** Restrict `clusters_to_use` to populations relevant to the question and provide a supported root when available. Changing the root changes the pseudotime reference, not the expression data.
- **Graph versus genes:** `num_dim`, `resolution`, and `k` affect the constructed representation and graph; `q_value_threshold` filters graph-associated genes after learning. A different gene cutoff does not repair an inappropriate trajectory.

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
    -t trajectory
```

## Outputs

Files are written under `RESULTS_DIR/monocle3/`.

| Files | Contents |
|---|---|
| `cds_object.rds` | Monocle3 cell_data_set with graph and ordering |
| `trajectory_summary.rds` | Run summary, including manual/automatic root mode |
| `trajectory_dependent_genes.csv` | Graph-test gene statistics |
| `plot/` | Trajectory, pseudotime, and available gene visualizations |

## Inspect the results

[![Principal graphs overlaid on the cluster-colored Monocle3 embedding.](../assets/examples/trajectory.png)](../assets/examples/trajectory.png)

**Example report.** This example contains multiple disconnected components. Do not interpret them as one continuous developmental path. Review the selected cell population, root, and finite pseudotime coverage before drawing lineage conclusions.


!!! note "Note"
    - **Inspect:** Check root-selection mode in `trajectory_summary.rds`, then review the trajectory and pseudotime views in `plot/`. Examine graph connectivity, represented populations, and cells without finite pseudotime.
    - **Decide:** Resolve an inappropriate cell subset or root before interpreting progression. Use `trajectory_dependent_genes.csv` for retained graph-associated genes, not as a replacement for a group-level differential-expression test.
    - **Interpretation:** Pseudotime represents an inferred ordering of transcriptional states,
      not chronological time. The selected root determines the direction of that ordering.
    - **Assessment:** Evaluate the biological relevance of the selected cells, the root recorded
      in the run summary, and the connectivity of the learned graph before interpreting trajectories.
    - **Limitations:** Disconnected partitions may contain cells with infinite pseudotime, and
      an unsuccessful graph test may produce an empty gene table. Inspect the stage log and
      ordering results before drawing conclusions about progression or associated genes.
