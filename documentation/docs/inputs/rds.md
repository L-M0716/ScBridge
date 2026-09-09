# Seurat RDS Input

RDS mode resumes analysis from an existing Seurat object at a supported processing stage.

## Input requirements

- **Input mode:** `-I rds`.
- **Input path (`-P`):** A Seurat RDS file or a supported filtering-stage directory.
- **Entry stage (`-G`):** The completed processing stage: `filtering`, `normalization`, `clustering`, or `annotation`.
- **Configuration (`-C`):** A YAML file containing organism and analysis settings.
- **Metadata (`-S`):** A tab-separated sample metadata table.
- **Annotation markers (`-M`):** A marker table appropriate for the organism and tissue.
- **Output (`-R`):** A writable results directory.
- **Analysis target (`-t`):** The requested endpoint or module; use an explicit target compatible with the entry stage.
- **Object contents:** The assays, reductions, and cell metadata required by the selected subsequent steps.
- **Sample identifiers:** For multi-sample filtering input, directory names matching `sample_id` in the metadata.

## Supported entry stages

<div class="rds-stages-table" markdown="1">

| Input object | Entry stage (`-G`) | Subsequent processing |
|---|---|---|
| Filtered Seurat object | `filtering` | Normalization; Batch Correction and Feature Selection; Clustering and Cell Annotation; Downstream Analysis |
| Normalized Seurat object | `normalization` | Batch Correction and Feature Selection; Clustering and Cell Annotation; Downstream Analysis |
| Clustered Seurat object | `clustering` | Cell Annotation; Downstream Analysis |
| Annotated Seurat object | `annotation` | Downstream Analysis |

</div>

The selected target determines which subsequent stages run. Downstream analysis includes differential
expression, GO/KEGG enrichment, Monocle3 trajectory analysis, and CellChat.

## RDS input

`-G, --input-stage` describes the processing already completed; `-P, --stage-path` identifies the
input object. Snakemake resolves the required steps between that entry stage and the requested target.

<div class="rds-paths-table" markdown="1">

| Entry stage | Input path (`-P`) |
|---|---|
| `filtering` | One RDS file for a single sample, or a directory containing `<sample_id>/seurat_filtered.rds` for each sample |
| `normalization` | A normalized Seurat RDS file |
| `clustering` | A clustered Seurat RDS file |
| `annotation` | An annotated Seurat RDS file |

</div>

!!! note "Stage compatibility"
    Choose the stage that matches the object's contents, not just its filename. Earlier analysis
    targets cannot be recovered by selecting a later-stage object. For multi-sample filtering input,
    directory names must match the metadata `sample_id` values.

## Processing flow

```text
Existing Seurat RDS
  ↓
Completed input stage (-G)
  ↓
Required subsequent steps
  ↓
Requested target (-t)
```

## Example command

Run from the repository root with `scRNA_seq.sif` available there. This example starts from a
clustered object and runs cell annotation. Use `differential`, `enrichment`, `trajectory`, or
`cellchat` for the corresponding downstream analysis; see [Workflow Targets](../reference/targets.md).

- **Paths:** Use container-accessible paths; bind directories outside the repository separately.
- **Preview:** Append `-n` for a dry run. Setup may still generate runtime configuration and standardization files.

```bash
apptainer exec \
  -B "$PWD":/opt/scRNA_workflow \
  scRNA_seq.sif \
  bash /opt/scRNA_workflow/run_workflow \
    -I rds \
    -G clustering \
    -P /opt/scRNA_workflow/rds/seurat_clustered.rds \
    -C /opt/scRNA_workflow/config/config.yaml \
    -S /opt/scRNA_workflow/config/samples.tsv \
    -M /opt/scRNA_workflow/config/markers.tsv \
    -R /opt/scRNA_workflow/results_rds \
    -t annotation
```
