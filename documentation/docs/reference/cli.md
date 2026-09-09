<div class="reference-tables" markdown="1">

# Command Line Reference

## General syntax

```text
bash /opt/scRNA_workflow/run_workflow -I fastq|matrix|rds [mode-specific input] \
  -C CONFIG -S METADATA -M MARKERS -R RESULTS [-t TASK]
```

## Arguments

| Argument | Scope | Purpose |
|---|---|---|
| `-I, --input-mode` | All | Select `fastq`, `matrix`, or `rds` |
| `-F, --fastq-dir` | FASTQ | FASTQ directory |
| `-D, --matrix-dir` | Matrix | Matrix directory |
| `-G, --input-stage` | RDS | `filtering`, `normalization`, `clustering`, or `annotation` |
| `-P, --stage-path` | RDS | RDS file; a per-sample directory is supported only for the `filtering` entry stage |
| `-C, --config` | All | Workflow configuration YAML |
| `-S, --metadata` | All | Metadata TSV |
| `-M, --markerlist` | All | Marker gene table |
| `-R, --results` | All | Output directory |
| `--star-index` | FASTQ | STAR index directory |
| `-t, --task` | All | Select workflow target |
| `-n, --dry-run` | All | Preview planned jobs |
| `--use-conda` | All | Use the Conda environment specified by each rule; enabled automatically by the launcher |
| `-h, --help` | All | Print launcher help |

## Targets and aliases

See [Workflow Targets](targets.md) for accepted targets and per-stage compatibility.
The default task is `all`; use an explicit target for RDS runs.

Accepted aliases include `seurat` for `rds`, `--samples` for `--metadata`,
`filtered` for `filtering`, `normalized` for `normalization`, `clustered` for
`clustering`, `annotated` for `annotation`, `de` for `differential`, and `enrich` for `enrichment`.

## Container execution

Run the Bash launcher inside the container as shown in [Quick Start](../quickstart.md).
`$PWD` identifies the host directory being mounted; analysis arguments use container-side paths.
`-n` previews jobs but does not suppress all setup file writes.

</div>
