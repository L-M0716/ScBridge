<div class="reference-tables" markdown="1">

# Troubleshooting

Start with the first failing job and its stage log. Later missing-output errors may be consequences
of that earlier failure.

## Container startup

| Symptom | Check |
|---|---|
| Image not found | SIF path relative to the current host directory |
| Apptainer command unavailable | Host runtime installation or module loading |
| Mount failure or permission denied | Host container permissions and fakeroot support |
| File missing inside container | Bind mounts and container-side paths |

Test startup with `apptainer exec scRNA_seq.sif echo ok`. Resolve startup failures before testing
analysis parameters. A mount error alone does not establish that the image is damaged.

## Input paths and identifiers

| Symptom | Check |
|---|---|
| FASTQ not found | `<run_id>_1.fastq.gz` and `<run_id>_2.fastq.gz`; metadata run IDs |
| Invalid matrix layout | `<sample_id>/filtered_feature_bc_matrix/` and the three gzip files |
| Filtering RDS missing | `<sample_id>/seurat_filtered.rds`; one file only for a single sample |
| Later-stage RDS error | Supply an actual RDS file, not a directory, for normalization/clustering/annotation |
| Reference failure | Matching FASTA, GTF, and index; files visible inside the container |

## Target selection

Select a target listed in [Workflow Targets](targets.md) and ensure that it is compatible
with the input stage. For RDS input, specify the required analysis target explicitly.

## QC and annotation

| Symptom | Check |
|---|---|
| Excessive cell loss | Per-sample QC distributions, count/gene thresholds, mitochondrial pattern |
| No useful subtype assignments | Species filter, matched marker genes, subtype-parent mapping |
| Many unassigned cell-type labels | SingleR scores, reference coverage, cluster heterogeneity |
| Unexpected CellChat labels | Contents of `cell_type_column`, including externally supplied RDS metadata |

## Empty downstream results

| Analysis | Check |
|---|---|
| Differential expression | Two samples per group, cell counts per cell type, control label, counts layer |
| Enrichment | DEG selection, mapped Entrez count, offline KEGG availability, summary status |
| Monocle3 | Root label match, selected clusters, graph connectivity, graph-test log |
| CellChat | At least two cell types, cell-group sizes, skipped small groups, status object contents |

!!! note "Note"
    A skipped analysis does not establish a negative biological result. Consult the corresponding
    analysis page for input requirements and guidance on interpreting its outputs.

</div>
