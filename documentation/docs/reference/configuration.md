<div class="reference-tables" markdown="1">

# Configuration Reference

Use `-C` to supply the YAML file. Tables list example parameter values; study-specific values
and paths should match the dataset and execution environment.

## File structure

Analysis settings are grouped under `qc`, `normalization`, `integration`, `clustering`,
`differential`, and `downstream`. FASTQ-specific settings use `reference` and `alignment`.

The launcher writes a runtime configuration that replaces input and output paths from the command-line
arguments. It also supplies input mode and RDS stage information. Paths passed to the launcher must
be accessible inside the container.

## Global

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `organism` | `human` | Species for annotation, enrichment, and CellChat: human or mouse. |
| `platform` | `10x_v3` | FASTQ library configuration: 10x_v3 or BD_v1. |

</div>

## Quality control

Keys in this table are within `qc`.

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `min_genes` | `200` | Minimum detected genes per cell, inclusive. |
| `max_genes` | `5000` | Maximum detected genes per cell, inclusive. |
| `min_counts` | `700` | Minimum total counts per cell, inclusive. |
| `max_mito` | `10` | Maximum mitochondrial percentage, inclusive. |
| `mito_pattern` | `^MT-` | Case-sensitive pattern for mitochondrial gene names. |
| `doublet_method` | `scDblFinder` | Recorded method label; the current script always calls scDblFinder. |
| `expected_doublet_rate` | `0.06` | Expected doublet fraction passed as dbr; 0.06 means 6%. |
| `gene_filtering.min_cells` | `7` | Minimum cells expressing a retained gene; also passed to SCTransform. |
| <code>gene_<wbr>filtering.<wbr>exclude_<wbr>genes</code> | <code>["^Rp[sl]","^Mt-","Malat1"]</code> | Case-insensitive gene-name exclusion patterns. |

</div>

## Normalization

Keys in this table are within `normalization`.

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `n_variable_genes` | `3000` | Consensus feature count during sample aggregation, not per-sample SCT feature count. |
| `vars_to_regress` | `["percent.mt"]` | Metadata variables regressed during SCTransform. |
| `seed` | `42` | Shared random seed for several core stages. |

</div>

## Integration

Keys in this table are within `integration`.

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `method` | `seurat` | seurat, harmony, or combat. |
| `batch_key` | `batch` | Metadata column defining batches. |
| `n_pcs` | `30` | Dimensions for integration, neighbor graph, and UMAP. |
| `n_features` | `3000` | Fallback feature count when variable features are absent. |

</div>

## Clustering

Keys in this table are within `clustering`.

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `resolution.default` | `0.8` | Graph-clustering resolution. |
| `n_neighbors` | `20` | FindNeighbors k.param. |
| `algorithm` | `Leiden` | Leiden or Louvain; unknown values currently fall back to Louvain. |
| `umap.min_dist` | `0.3` | UMAP min.dist. |
| `umap.spread` | `1` | UMAP spread. |

</div>

## Differential expression

Keys in this table are within `differential`.

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `group_by_column` | `group` | Biological-group column in sample and cell metadata. |
| `control_group` | `Adult peripheral blood` | Reference-group label matching the sample metadata. |
| `logfc_threshold` | `0.25` | DESeq2 lfcThreshold and annotation marker-report threshold. |
| `min_pct` | `0.1` | Annotation marker-report detection fraction; not used in the DESeq2 test. |
| `min_cells_group` | `20` | Minimum cells per group within each cell type. |
| `visualization.top_n` | `20` | Number of genes in selected summaries. |
| <code>visualization.<wbr>adj_<wbr>pval_<wbr>cutoff</code> | `0.05` | Adjusted-P cutoff for visualization. |
| <code>visualization.<wbr>logfc_<wbr>threshold</code> | `0.5` | Effect-size cutoff for visualization. |

</div>

## Trajectory

Keys in this table are within `downstream.trajectory`.

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `clusters_to_use` | `all` | all or comma-separated Seurat cluster IDs. |
| `root_cluster` | `""` | Optional Seurat cluster used when no cell-type reference is applied. |
| `root_cell_type` | `""` | Optional annotation label; takes precedence when cell-type metadata are available. |
| `num_dim` | `50` | Monocle3 preprocessing dimensions. |
| `resolution` | `0.001` | Monocle3 clustering resolution. |
| `k` | `20` | Monocle3 neighbor parameter. |
| `q_value_threshold` | `0.05` | Graph-test gene selection cutoff. |

</div>

`""` denotes an empty string: no root cell type or cluster is specified. When both root settings
are empty, the workflow selects a root automatically. See
[Trajectory root specification](../downstream/monocle3.md#root-selection) for the selection rules.

## Communication

Keys in this table are within `downstream.cellchat`.

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `db` | `all` | Full database (all) or a supported subset. |
| `min_cells` | `10` | Minimum cell-group size for communication filtering. |
| `cell_type_column` | `singler_labels` | Metadata column used for communication identities. |
| `pathways_of_interest` | `[]` | Optional pathway names for visualization. |
| `top_n_pathways` | `20` | Number of pathways displayed. |
| `n_patterns_outgoing` | `4` | Requested outgoing signaling patterns. |
| `n_patterns_incoming` | `4` | Requested incoming signaling patterns. |

</div>

`[]` denotes an empty list: no specific pathways are requested for visualization. The report
then uses the first `top_n_pathways` entries in the inferred pathway list.

## Enrichment

Keys in this table are within `downstream.pathway_enrichment`.

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `databases.mouse_db_path` | <code>Species-specific RDS path</code> | Species-specific offline KEGG RDS path. |
| `databases.human_db_path` | <code>Species-specific RDS path</code> | Species-specific offline KEGG RDS path. |
| `parameters.pval_cutoff` | `0.05` | Input DEG adjusted-P threshold and enrichment P-value cutoff. |
| `parameters.qval_cutoff` | `0.1` | Enrichment q-value cutoff. |
| <code>parameters.<wbr>logfc_<wbr>threshold</code> | `0.5` | Absolute DEG log2 fold-change filter. |
| `parameters.min_gene_size` | `10` | Minimum enrichment term size. |
| `parameters.max_gene_size` | `500` | Maximum enrichment term size. |
| <code>visualization.<wbr>top_<wbr>n_<wbr>terms</code> | `20` | Terms displayed per plot. |
| `visualization.plot_width` | `12` | Plot width in inches. |
| <code>visualization.<wbr>plot_<wbr>height</code> | `8` | Plot height in inches. |

</div>

## FASTQ alignment

Keys in this table are within `alignment`.

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `read_length` | `90` | Reference-index read length; sjdbOverhang uses read length minus one. |
| `soloBarcodeReadLength` | `28` | 10x v3 barcode-read length. |
| `soloCBlen` | `16` | 10x v3 cell-barcode length. |
| `soloUMIlen` | `12` | 10x v3 UMI length. |
| `bd_v1.soloCBposition` | <code>0_<wbr>0_<wbr>0_<wbr>8 0_<wbr>21_<wbr>0_<wbr>29 0_<wbr>43_<wbr>0_<wbr>51</code> | BD v1 barcode segments for CB_UMI_Complex. |
| `bd_v1.soloUMIposition` | `0_52_0_59` | BD v1 UMI position. |
| `bd_v1.soloStrand` | `Forward` | BD v1 transcript strand setting. |

</div>

## Paths and resources

| Setting | Role | Runtime behavior |
|---|---|---|
| `input.fastq_dir` | FASTQ input directory | Replaced by `-F` |
| `input.matrix_dir` | Matrix input directory | Supplied by `-D` |
| `input.samples_table` | Sample metadata | Replaced by `-S` |
| `input.cell_cycle_genes` | Cell-cycle gene RDS | Species-matched resource; not switched by organism |
| <code>reference.<wbr>paths.<wbr>genome_<wbr>fa</code> | Genome FASTA | Required for FASTQ reference processing |
| `reference.paths.gtf` | Genome annotation | Same assembly as the FASTA |
| <code>reference.<wbr>paths.<wbr>star_<wbr>index</code> | STAR index directory | Replaced by `--star-index`, or `RESULTS_DIR/STAR_index` |
| <code>reference.<wbr>paths.<wbr>whitelists</code> | Library barcode resources | 10x v3 file or three BD v1 files |
| <code>database_<wbr>resources.<wbr>reference_<wbr>databases.<wbr>markerlist</code> | Annotation marker table | Replaced by `-M` |
| `output.results` | Stage output directories | Rebased under `-R` |
| `script_paths` | R script locations | Must match the mounted source tree |
| `conda.doublet_detect` | Environment prefix | Must exist inside the execution environment |

The genome FASTA and GTF must correspond to the same assembly. Any reused STAR index must be
compatible with these reference files.

## Computing resources

Stage-specific thread and memory settings are grouped under `resources` and read by the
Snakemake rules. The examples below are resource declarations, not measured usage or minimum
hardware requirements. Actual requirements depend on dataset size and the analysis method.

<div class="parameter-table" markdown="1">

| Parameter | Example | Description |
|---|---|---|
| `default_threads` | `8` | Requested rule threads |
| `default_mem_gb` | `32` | Memory in GB |
| `alignment.star.threads` | `18` | Requested rule threads |
| `alignment.star.mem_gb` | `128` | Memory in GB |
| <code>alignment.<wbr>star.<wbr>index_<wbr>threads</code> | `8` | Requested rule threads |
| <code>normalization.<wbr>sctransform.<wbr>threads</code> | `8` | Requested rule threads |
| <code>normalization.<wbr>sctransform.<wbr>mem_<wbr>gb</code> | `32` | Memory in GB |
| <code>batch_<wbr>correction.<wbr>correction.<wbr>threads</code> | `12` | Requested rule threads |
| <code>batch_<wbr>correction.<wbr>correction.<wbr>mem_<wbr>gb</code> | `128` | Memory in GB |
| `clustering.pca.threads` | `8` | Requested rule threads |
| `clustering.pca.mem_gb` | `32` | Memory in GB |
| `annotation.auto.threads` | `8` | Requested rule threads |
| `annotation.auto.mem_gb` | `32` | Memory in GB |
| `differential.run.threads` | `8` | Requested rule threads |
| `differential.run.mem_gb` | `32` | Memory in GB |
| `differential.viz.threads` | `4` | Requested rule threads |
| `differential.viz.mem_gb` | `16` | Memory in GB |
| <code>downstream.<wbr>monocle3.<wbr>threads</code> | `8` | Requested rule threads |
| <code>downstream.<wbr>monocle3.<wbr>mem_<wbr>gb</code> | `32` | Memory in GB |
| <code>downstream.<wbr>cellchat.<wbr>threads</code> | `8` | Requested rule threads |
| <code>downstream.<wbr>cellchat.<wbr>mem_<wbr>gb</code> | `32` | Memory in GB |
| <code>downstream.<wbr>pathway_<wbr>enrichment.<wbr>threads</code> | `4` | Requested rule threads |
| <code>downstream.<wbr>pathway_<wbr>enrichment.<wbr>mem_<wbr>gb</code> | `16` | Memory in GB |
| <code>doublet_<wbr>detection.<wbr>threads</code> | `8` | Requested rule threads |
| `doublet_detection.mem_gb` | `32` | Memory in GB |
| `filtering.threads` | `8` | Requested rule threads |
| `filtering.mem_gb` | `32` | Memory in GB |

</div>

The launcher currently passes `--cores all` and `--use-conda` to Snakemake. Individual R scripts
may use fewer threads than their rule allocation.

<span id="fixed-implementation-settings"></span>

!!! note "Note"
    Only settings read by the workflow rules or analysis scripts affect execution. Adding a YAML
    key does not introduce an alternative method or expose a script-level constant.
    Method-specific settings are described in [Normalization](../workflow/normalization.md),
    [Cell Annotation](../workflow/annotation.md), [Differential Expression](../downstream/differential.md),
    and [Cell-Cell Communication](../downstream/cellchat.md).

</div>
