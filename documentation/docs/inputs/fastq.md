# FASTQ Input

FASTQ mode starts from paired-end sequencing reads, performs read-quality assessment, and generates
gene-count matrices with STARsolo before entering the shared ScBridge workflow.

## Input requirements

- **Input mode:** `-I fastq`.
- **Input path (`-F`):** A directory of gzip-compressed paired-end FASTQ files named by run ID.
- **Configuration (`-C`):** A YAML file containing organism and analysis settings.
- **Metadata (`-S`):** A tab-separated sample metadata table.
- **Annotation markers (`-M`):** A marker table appropriate for the organism and tissue.
- **Output (`-R`):** A writable results directory.
- **Analysis target (`-t`):** The requested endpoint or module; `all` runs the full workflow from the selected input stage.
- **Read metadata:** `sample_id`, `batch`, `group`, and `run_ids` columns, in that order; comma-separated run IDs for samples with multiple runs.
- **Reference and library settings:** Matching genome FASTA and GTF paths, a supported library platform, and barcode resources in the configuration. `--star-index` optionally specifies an existing compatible index directory.

### Library configuration

The current alignment rule supports these exact `platform` values:

| Platform | STARsolo mode | Barcode and UMI settings |
|---|---|---|
| `10x_v3` | `CB_UMI_Simple` | Uses `alignment.soloCBlen`, `soloUMIlen`, and `soloBarcodeReadLength`; the supplied defaults are 16, 12, and 28 |
| `BD_v1` | `CB_UMI_Complex` | Uses `alignment.bd_v1.soloCBposition`, `soloUMIposition`, and `soloStrand` |

For 10x v3, `reference.paths.whitelists.10xv3` selects the whitelist; the bundled example is
`3M-february-2018_TRU.txt`. For BD v1, `reference.paths.whitelists.bd_v1` supplies the three
`BD_CLS1_V1.txt`, `BD_CLS2_V1.txt`, and `BD_CLS3_V1.txt` files.
Both branches pass Read 2 as cDNA and Read 1 as barcode/UMI input.
Confirm that the read structure and resources match your library; these settings do not imply support
for every 10x or BD chemistry.

## FASTQ layout

The directory supplied with `-F` contains one read pair for each run ID in the metadata:

```text
fastq_files/
  RUN001_1.fastq.gz
  RUN001_2.fastq.gz
  RUN002_1.fastq.gz
  RUN002_2.fastq.gz
```

Here, `RUN001` and `RUN002` are illustrative identifiers. Replace them with the run IDs in your
metadata. The current rule expects the `_1.fastq.gz` and `_2.fastq.gz` suffixes.

## FastQC

ScBridge runs FastQC on both reads of each run before STARsolo. HTML reports and ZIP archives are
written to `RESULTS_DIR/read_data/runs/RUN_ID/`.
FastQC reports read-quality metrics; this step does not trim reads or replace cell-level quality control.

## STAR and STARsolo

ScBridge uses STARsolo for single-cell FASTQ alignment and gene-count matrix generation.
STARsolo is integrated into the STAR executable, so no separate STARsolo installation is required.

| Operation | Key STAR/STARsolo option |
|---|---|
| Reference index construction | `--runMode genomeGenerate` |
| 10x v3 processing | `--soloType CB_UMI_Simple` |
| BD v1 processing | `--soloType CB_UMI_Complex` |

!!! note "STARsolo settings"
    ScBridge supplies these options to `STAR`; they are not standalone commands.
    Barcode, UMI, and whitelist settings follow the library configuration. Both library modes use
    `GeneFull` counting and `EmptyDrops_CR` cell calling, as defined in the alignment rule.

For algorithm and parameter details, see the
[official STARsolo documentation](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md).

### Reference index behavior

The workflow includes a `build_star_index` rule that runs `STAR --runMode genomeGenerate` using
`reference.paths.genome_fa` and `reference.paths.gtf`. STARsolo depends on its `Genome` and `SA`
outputs, allowing Snakemake to schedule index construction when those files are missing.

Use `--star-index` to select the index directory. If omitted, the launcher uses
`RESULTS_DIR/STAR_index`. Reference files must be accessible inside the container, and the destination
must be writable when building a new index.

!!! note "Reusing an index"
    The build rule skips construction when both `Genome` and `SA` exist. This file-existence check
    does not validate reference compatibility or repair an incomplete index. Reuse a complete index
    matching your genome and annotation; select a new, empty directory when changing references.

## Processing flow

```text
FASTQ
  ↓
FastQC
  ↓
STARsolo
  ↓
Gene × cell count matrix
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

Run from the ScBridge repository root with the SIF image available there. This example uses the
same input layout as [Quick Start](../quickstart.md); set the reference paths in
`config/config.yaml` to match your files.

- **Paths:** Use container-accessible paths; bind directories outside the repository separately.
- **Preview:** Append `-n` for a dry run. Setup may still generate runtime configuration and standardization files.

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
    -t all
```
