# Pipeline Overview

ScBridge resolves analysis dependencies from the input mode and the requested target.
RDS entry stages bypass completed processing; they do not restart at raw-read or cell-level QC.

## Input entry points

| Input mode | Starting data | Entry point |
|---|---|---|
| FASTQ | Paired-end reads | FastQC and STARsolo, then count matrices |
| Matrix | 10x-style count matrices | Cell-level QC |
| RDS | Supported Seurat object | After filtering, normalization, clustering, or annotation |

See [Input Data](../inputs/index.md) for file requirements.

## Core workflow

The sequence below lists analysis operations from cell-level QC to cell-type annotation.
Count matrices are the input to this sequence; FASTQ preprocessing is described above.

```text
Quality-control assessment
  ↓
Doublet detection
  ↓
Cell and gene filtering
  ↓
Normalization and consensus-feature selection
  ↓
Batch correction and integration
  ↓
Variable-feature assessment
  ↓
Neighbor-graph construction and clustering
  ↓
UMAP visualization
  ↓
Reference-based cell-type annotation
  ↓
Marker-based subtype refinement
```

## RDS continuation

<div class="reference-tables" markdown="1">

| Completed stage | Subsequent processing |
|---|---|
| `filtering` | Normalization, Integration, Feature Selection, Clustering, Annotation |
| `normalization` | Integration, Feature Selection, Clustering, Annotation |
| `clustering` | Annotation |
| `annotation` | Downstream Analyses |

</div>

Use `-G` for the completed input stage and `-t` for the requested endpoint.
The table summarizes the continuation sequence; the selected target determines where processing
stops. See [RDS target restrictions](../reference/targets.md#rds-target-restrictions) for accepted
target names.

## Downstream branches

```text
Annotated object
  ├── Differential expression (DESeq2 pseudobulk)
  │     └── Functional enrichment (GO / KEGG)
  ├── Trajectory inference (Monocle3)
  └── Cell-cell communication (CellChat)
```

Trajectory and communication inference do not depend on enrichment.
Enrichment requires differential-expression results.

## Stage-aware execution

`-t` selects a workflow target, not an isolated script. Snakemake schedules missing dependencies
and can reuse complete, up-to-date outputs. Input stage, configuration, and target all affect the planned jobs.

See [Workflow Targets](../reference/targets.md) for accepted targets and RDS-stage compatibility. Append `-n` to preview jobs; setup may still generate runtime files.
