<p align="center">
  <img src="ScBridge_logo.png" alt="ScBridge - Bridging single-cell data to insights" width="750">
</p>

<p align="center">
  <strong>A containerized Snakemake workflow for end-to-end single-cell RNA-seq analysis.</strong>
</p>

ScBridge supports FASTQ, 10x matrix, and Seurat RDS inputs, combining
hierarchical cell-type annotation with comprehensive downstream analysis.
## Workflow Overview

![ScBridge workflow](ScBridge_workflow.png)
## Highlights

<ul>
  <li>Three entry points: FASTQ, 10x matrix, and Seurat RDS.</li>
  <li>Stage-aware RDS input from filtering, normalization, clustering, or annotation.</li>
  <li>FastQC and STARsolo processing for raw FASTQ data.</li>
  <li>SingleR/celldex major-lineage annotation followed by marker-based subtype refinement.</li>
  <li>Differential expression, enrichment, Monocle3, and CellChat analyses.</li>
  <li>Reproducible dependency management through one Apptainer SIF image.</li>
  <li>Snakemake schedules only the jobs required for the selected target and reuses complete, up-to-date outputs.</li>
</ul>

## Project at a Glance

| Component | Implementation |
|---|---|
| Workflow engine | Snakemake |
| Containerization | Apptainer with a SIF image |
| FASTQ processing | FastQC and STARsolo |
| Core analysis | Seurat-based QC, filtering, normalization, integration, and clustering |
| Annotation | SingleR/celldex major types plus marker-based subtype refinement |
| Downstream analysis | Differential expression, GO/KEGG enrichment, Monocle3, and CellChat |

## Installation

ScBridge requires a Linux system with Apptainer. R, Python, Conda, and workflow
packages are provided inside the image and do not need to be installed on the
host.

Build the image from the supplied definition file:

```bash
cd <PROJECT_DIR>

apptainer build --fakeroot scRNA_seq.sif \
  build/scRNA_seq.def
```

## Quick Start

The general command is:

```text
run_workflow -I fastq|matrix|rds [mode-specific input] \
  -C CONFIG -S METADATA -M MARKERS -R RESULTS [-t TASK]
```

### FASTQ Input

```bash
apptainer exec \
  -B <PROJECT_DIR>:/opt/scRNA_workflow \
  scRNA_seq.sif \
  bash /opt/scRNA_workflow/run_workflow \
    -I fastq \
    -F <FASTQ_DIR> \
    -C <CONFIG_FILE> \
    -S <METADATA_FILE> \
    -M <MARKER_FILE> \
    -R <RESULTS_DIR> \
    --star-index <STAR_INDEX_DIR> \
    -t all
```

If the STAR index does not exist, ScBridge can build it from the genome FASTA
and GTF paths defined in `config.yaml`.

### Matrix Input

```bash
apptainer exec \
  -B <PROJECT_DIR>:/opt/scRNA_workflow \
  scRNA_seq.sif \
  bash /opt/scRNA_workflow/run_workflow \
    -I matrix \
    -D <MATRIX_DIR> \
    -C <CONFIG_FILE> \
    -S <METADATA_FILE> \
    -M <MARKER_FILE> \
    -R <RESULTS_DIR> \
    -t all
```

### RDS Input

```bash
apptainer exec \
  -B <PROJECT_DIR>:/opt/scRNA_workflow \
  scRNA_seq.sif \
  bash /opt/scRNA_workflow/run_workflow \
    -I rds \
    -G <INPUT_STAGE> \
    -P <RDS_PATH> \
    -C <CONFIG_FILE> \
    -S <METADATA_FILE> \
    -M <MARKER_FILE> \
    -R <RESULTS_DIR> \
    -t all
```

For RDS mode, `-G` accepts `filtering`, `normalization`, `clustering`, or
`annotation`. The `-P` argument may point to a supported RDS file or a
multi-sample stage directory.

Replace each uppercase placeholder with an absolute path visible inside the
container. With the mount shown above, files under `<PROJECT_DIR>` are available
under `/opt/scRNA_workflow`.

## Inputs

| Argument | Purpose |
|---|---|
| `-I, --input-mode` | Select `fastq`, `matrix`, or `rds` input |
| `-F, --fastq-dir` | FASTQ directory for FASTQ mode |
| `-D, --matrix-dir` | Directory containing one 10x matrix folder per sample |
| `-G, --input-stage` | Stage represented by an RDS input |
| `-P, --stage-path` | RDS file or stage directory |
| `-C, --config` | Workflow configuration YAML |
| `-S, --metadata` | Tab-separated sample metadata |
| `-M, --markerlist` | Species-aware marker gene table |
| `-R, --results` | Output directory for the analysis |

Sample identifiers in the metadata must match the FASTQ, matrix, or RDS sample
names used by the selected input mode.

## Tasks

Select a target with `-t, --task`:

| Task | Behaviour |
|---|---|
| `qc` | Stop after quality control |
| `filtering` | Stop after filtering |
| `normalization` | Stop after normalization |
| `clustering` | Stop after clustering |
| `annotation` | Stop after cell-type annotation |
| `differential` | Run differential expression analysis |
| `enrichment` | Run functional enrichment analysis |
| `trajectory` | Run Monocle3 trajectory analysis |
| `cellchat` | Run CellChat communication analysis |
| `downstream` | Run differential, enrichment, trajectory, and CellChat modules |
| `all` | Run from the selected input point to all available final outputs |

## Main Outputs

Each analysis is isolated under the directory supplied with `-R`:

```text
<RESULTS_DIR>/
├── qc/
├── doublets/
├── filtering/
├── normalization/
├── batch_correction/
├── feature_selection/
├── clustered/
├── annotation/
│   ├── seurat_final_annotated.rds
│   ├── reports/
│   └── spreadsheets/
├── differential/
│   ├── reports/
│   └── spreadsheets/
├── enrichment/
│   ├── enrichment_summary.csv
│   ├── enrichment_results_all.rds
│   └── visualizations
├── monocle3/
│   ├── cds_object.rds
│   ├── trajectory_summary.rds
│   ├── trajectory_dependent_genes.csv
│   └── plot/
├── cellchat/
│   ├── plots/
│   └── spreadsheets/
├── reports/
│   └── standardization_report.txt
└── runtime/
    └── input manifests used by the workflow
```

Earlier-stage directories are created only when those stages are executed in
the current analysis.

## Documentation

The complete guide covers input preparation, configuration, command-line use,
pipeline stages, output interpretation, and troubleshooting:

* [Online documentation](https://scbridge.readthedocs.io/)
* [Documentation source](docs/index.md)
* [Open the local documentation](OPEN_SCBridge_DOCUMENTATION.html)

## Citation And License

Citation information will be added with the first public release. See
[`LICENSE`](LICENSE) for the project licence.
