<div class="quickstart-guide" markdown="1">

# Quick Start

## Choose an input mode

| Mode | Starting point | Required input arguments |
|---|---|---|
| `fastq` | Starting from raw FASTQ reads | `-F FASTQ_DIR` |
| `matrix` | Starting from existing 10x-style count matrices | `-D MATRIX_DIR` |
| `rds` | Continuing from supported Seurat objects | `-G STAGE -P PATH` |

## Prepare the input files

Prepare only the inputs for your chosen mode. These layouts are examples, not required directory names.
Use paths matching your files in the command-line arguments and configuration.
The examples below assume that the input directories and a `config/` directory are in the ScBridge
repository root, alongside `run_workflow` and `scRNA_seq.sif`.

=== "FASTQ"

    ```text
    config/
      config.yaml
      samples.tsv
      markers.tsv
    fastq_files/
      <run_id>_1.fastq.gz
      <run_id>_2.fastq.gz
    ```

    Set the genome FASTA, GTF, and library-specific barcode resources in the configuration.
    See [FASTQ Input](inputs/fastq.md) for reference and library requirements.

=== "Matrix"

    ```text
    config/
      config.yaml
      samples.tsv
      markers.tsv
    matrix_files/
      <sample_id>/
        matrix.mtx.gz
        features.tsv.gz
        barcodes.tsv.gz
    ```

    Match the sample folders to the metadata table. See [Matrix Input](inputs/matrix.md)
    for supported matrix formats.

=== "RDS"

    ```text
    config/
      config.yaml
      samples.tsv
      markers.tsv
    filtering/
      <sample_id>/
        seurat_filtered.rds
    ```

    This example starts after filtering. The RDS structure depends on the selected stage:
    filtering, normalization, clustering, or annotation. See [RDS Input](inputs/rds.md)
    for the corresponding requirements.

## Preview the workflow

- **Preview:** `-n` performs a dry run without executing the analysis. ScBridge may still generate runtime configuration and standardization files during setup.
- **Run location:** Execute the command from the repository root. `$PWD` mounts this directory at `/opt/scRNA_workflow` inside the container.
- **Target:** All examples stop after annotation. The RDS example starts from filtering outputs; each mode uses a separate results directory.
- **Paths:** Adjust the examples to your files and use container-accessible paths for all workflow arguments.

=== "FASTQ"

    ```bash
    apptainer exec \
      -B "$PWD":/opt/scRNA_workflow \
      scRNA_seq.sif \
      bash /opt/scRNA_workflow/run_workflow \
        -I fastq \
        -F /opt/scRNA_workflow/fastq_files \
        -C /opt/scRNA_workflow/config/config.yaml \
        -S /opt/scRNA_workflow/config/samples.tsv \
        -M /opt/scRNA_workflow/config/markers.tsv \
        -R /opt/scRNA_workflow/results_fastq \
        -t annotation \
        -n
    ```

=== "Matrix"

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
        -t annotation \
        -n
    ```

=== "RDS"

    ```bash
    apptainer exec \
      -B "$PWD":/opt/scRNA_workflow \
      scRNA_seq.sif \
      bash /opt/scRNA_workflow/run_workflow \
        -I rds \
        -G filtering \
        -P /opt/scRNA_workflow/filtering \
        -C /opt/scRNA_workflow/config/config.yaml \
        -S /opt/scRNA_workflow/config/samples.tsv \
        -M /opt/scRNA_workflow/config/markers.tsv \
        -R /opt/scRNA_workflow/results_rds \
        -t annotation \
        -n
    ```

## Run the workflow

After confirming the expected jobs and paths, remove `-n` and run the same command.
The target `annotation` runs the required steps up to cell-type annotation.
Use `-t all` to include all downstream modules from the selected input point.
See [Command Line](reference/cli.md) for other analysis targets.

</div>
