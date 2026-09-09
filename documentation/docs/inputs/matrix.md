# 10x Matrix Input

Matrix mode accepts 10x-style count matrices for cell-level quality control and subsequent analysis.

## Input requirements

- **Input mode:** `-I matrix`.
- **Input path (`-D`):** A directory of per-sample 10x-style count matrices, barcode files, and feature files.
- **Configuration (`-C`):** A YAML file containing organism and analysis settings.
- **Metadata (`-S`):** A tab-separated sample metadata table.
- **Annotation markers (`-M`):** A marker table appropriate for the organism and tissue.
- **Output (`-R`):** A writable results directory.
- **Analysis target (`-t`):** The requested endpoint or module; `all` runs the full workflow from the selected input stage.
- **Sample identifiers:** Sample directory names matching `sample_id` in the metadata.
- **Matrix files:** `matrix.mtx.gz`, `barcodes.tsv.gz`, and `features.tsv.gz` within each sample's `filtered_feature_bc_matrix` directory.

## Directory structure

The directory supplied with `-D` must contain the following structure for each sample:

```text
matrix_files/
  sample01/
    filtered_feature_bc_matrix/
      matrix.mtx.gz
      barcodes.tsv.gz
      features.tsv.gz
  sample02/
    filtered_feature_bc_matrix/
      matrix.mtx.gz
      barcodes.tsv.gz
      features.tsv.gz
```

Replace `sample01` and `sample02` with the identifiers in your metadata.

## Matrix input

ScBridge reads the count matrix and its barcode and feature annotations directly, without FastQC
or STARsolo. Processing begins with cell-level quality control.

| Input file | Contents |
|---|---|
| `matrix.mtx.gz` | Sparse gene-by-cell counts |
| `barcodes.tsv.gz` | Cell barcode identifiers |
| `features.tsv.gz` | Gene or feature identifiers |

!!! note "Count data"
    Use unnormalized gene counts. ScBridge applies cell-level QC, doublet detection, and
    filtering to the supplied matrices, including matrices already filtered by the upstream software.

## Processing flow

```text
10x count matrix
  ↓
QC and filtering
  ↓
Normalization
  ↓
Batch correction and feature selection
  ↓
Clustering and annotation
  ↓
Requested downstream analysis
```

The selected target determines where processing stops. See [Pipeline Overview](../workflow/overview.md).

## Example command

Run from the repository root with `scRNA_seq.sif` available there. This example uses the directory
layout above and runs through the full analysis target.

- **Paths:** Use container-accessible paths; bind directories outside the repository separately.
- **Preview:** Append `-n` for a dry run. Setup may still generate runtime configuration and standardization files.

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
    -t all
```
