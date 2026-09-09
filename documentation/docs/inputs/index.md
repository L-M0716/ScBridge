# Input Data Overview

ScBridge supports three entry modes so analyses can begin from raw reads, existing count matrices,
or intermediate Seurat objects.

| Mode | Starting material | Typical use |
|---|---|---|
| `fastq` | Raw paired-end FASTQ files | End-to-end processing |
| `matrix` | 10x-style expression matrix | Re-analysis after count generation |
| `rds` | Seurat RDS file or stage directory | Resume from a supported workflow stage |

Current launcher behavior requires a configuration file, metadata table, marker table, and results directory
for all three modes.

## Human and mouse datasets

ScBridge supports human and mouse datasets. The `organism` setting in `config.yaml` (`human` or
`mouse`) determines the species-specific resources used by annotation, enrichment, and CellChat.
These resources are selected by the workflow; no separate database selection is required.

| Analysis resource | Human | Mouse |
|---|---|---|
| SingleR reference | `HumanPrimaryCellAtlasData` | `MouseRNAseqData` |
| GO annotation database | `org.Hs.eg.db` | `org.Mm.eg.db` |
| KEGG species code | `hsa` | `mmu` |
| CellChat database | `CellChatDB.human` | `CellChatDB.mouse` |

The organism is specified in the configuration, not inferred from the input data.

## Reference configuration

Only settings relevant to the requested analysis steps need review. Existing settings can be retained
when they already match the dataset; these are not additional selections required for every run.

| Analysis step | Dataset-specific resource | Configuration |
|---|---|---|
| FASTQ processing | Matching genome FASTA, GTF, and STAR index | `reference.paths`; `--star-index` |
| Quality control | Mitochondrial gene-name pattern | `qc.mito_pattern` |
| Normalization | Species-matched cell-cycle genes | `input.cell_cycle_genes` |
| Cell annotation | Species- and tissue-matched markers | `-M` |

!!! note "Species-specific inputs"
    The example configuration uses human resources. For mouse data, review the settings above for
    the steps you will run: changing `organism` does not update these paths or the mitochondrial
    pattern. Common mitochondrial symbol prefixes are `^MT-` (human) and `^mt-` (mouse);
    the pattern must match the gene names in your data.

## Marker tables

ScBridge accepts the standard columns `species`, `main_cell_type`, `sub_cell_type`, and `marker_gene`.
It also recognizes `official gene symbol` and `cell type`, mapping them to `marker_gene` and
`main_cell_type`, respectively.

| Marker species code | Organism |
|---|---|
| `Hs` | Human |
| `Mm` | Mouse |

Mixed-species marker tables are filtered by `organism`; shared entries such as `Mm Hs` are retained
for either species. Marker gene symbols must match the input data. This filtering does not perform
cross-species gene conversion. The bundled UCB marker table contains human entries.

**Subtype refinement requires subtype markers.** If `sub_cell_type` is absent, the script treats it
as empty and uses the table as a source of broad cell-type markers. To support refinement, provide
subtype labels, their parent cell types, and appropriate marker genes in the standard columns.
