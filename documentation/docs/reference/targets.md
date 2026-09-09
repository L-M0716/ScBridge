<div class="reference-tables" markdown="1">

# Workflow Targets

`-t` selects an output target and its dependencies. It does not mean that only one script runs.
There are no launcher targets for batch correction or feature selection alone.

## Analysis targets

| Target | Result | Immediate upstream requirement |
|---|---|---|
| `filtering` | Filtered cells and genes | QC metrics and doublet calls |
| `normalization` | Merged normalized object | Filtered sample objects |
| `clustering` | Clusters and UMAP | Normalization, integration, and feature stage |
| `annotation` | Final cell-type labels | Clustered object |
| `differential` | DESeq2 pseudobulk results | Annotated object and sample groups |
| `enrichment` | GO and KEGG results | Differential-expression table and annotated object |
| `trajectory` | Monocle3 graph and pseudotime | Annotated object |
| `cellchat` | Communication inference and reports | Annotated object |
| `all` | Request all processing after the selected entry point | Target mapping depends on the input mode and completed stage |

## RDS target restrictions

The table lists accepted values of `-t` for each RDS input stage. See
[RDS continuation](../workflow/overview.md#rds-continuation) for the processing sequence.

| Input stage | Explicit targets accepted |
|---|---|
| `filtering` | `normalization`, `clustering`, `annotation`, `differential`, `enrichment`, `trajectory`, `cellchat` |
| `normalization` | `clustering`, `annotation`, `differential`, `enrichment`, `trajectory`, `cellchat` |
| `clustering` | `annotation`, `differential`, `enrichment`, `trajectory`, `cellchat` |
| `annotation` | `differential`, `enrichment`, `trajectory`, `cellchat` |

The input object must contain the data needed by the selected analysis. Passing argument validation
does not validate every assay, reduction, or metadata field.

## Preview and reuse

Append `-n` to the complete command to inspect planned jobs without executing analysis.
Setup can copy RDS inputs and generate runtime files before the preview.

Use separate result directories for different input datasets or substantially different analyses.

</div>
