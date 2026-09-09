---
hide:
  - toc
---

<h1 class="scbridge-home-title">ScBridge: modular single-cell RNA-seq analysis</h1>

<div class="scbridge-hero">
  <img src="assets/scbridge_logo.png" alt="ScBridge logo">
</div>

ScBridge connects read processing, cell-type annotation, and downstream analysis in an Apptainer-packaged
Snakemake workflow. Start with FASTQ, 10x count matrices, or a supported Seurat RDS stage to obtain
processed objects, annotation evidence, and analysis reports.

## From input to analysis

<div class="pipeline-diagram" aria-label="ScBridge processing sequence and RDS entry points">
  <div class="pipeline-input">
    <span class="pipeline-label">INPUT</span>
    <div><a href="inputs/fastq.html">FASTQ</a><span>FastQC + STARsolo</span><span aria-hidden="true">&rarr;</span><strong>Count matrices</strong></div>
    <div><a href="inputs/matrix.html">10x Matrix</a><span aria-hidden="true">&rarr;</span><strong>Count matrices</strong></div>
    <div><a href="inputs/rds.html">Seurat RDS</a><span aria-hidden="true">&rarr;</span><strong>Completed stage</strong></div>
  </div>
  <section class="pipeline-phase"><div class="pipeline-phase-label">Stage 1<span>Preprocessing</span></div><div class="pipeline-phase-body">
  <ol class="pipeline-stages">
    <li><span class="pipeline-number">01</span><a href="workflow/qc_filtering.html">Quality control</a><p>Doublet detection<br>Cell filtering</p></li>
    <li><span class="pipeline-number">02</span><a href="workflow/normalization.html">Normalization</a><p>SCT normalization<br>Consensus features</p><a class="pipeline-resume" href="inputs/rds.html">Filtered RDS &uarr;</a></li>
    <li><span class="pipeline-number">03</span><a href="workflow/batch_correction.html">Integration</a><p>Batch correction<br><a href="workflow/feature_selection.html">Feature assessment</a></p><a class="pipeline-resume" href="inputs/rds.html">Normalized RDS &uarr;</a></li>
    <li><span class="pipeline-number">04</span><a href="workflow/clustering.html">Clustering</a><p>Graph-based clusters<br>UMAP visualization</p></li>
  </ol>
  </div></section>
  <section class="pipeline-phase"><div class="pipeline-phase-label">Stage 2<span>Annotation and analysis</span></div><div class="pipeline-phase-body">
  <ol class="pipeline-stages pipeline-annotation">
    <li><span class="pipeline-number">05</span><a href="workflow/annotation.html">Cell annotation</a><p>SingleR broad labels<br>UCell refinement</p><a class="pipeline-resume" href="inputs/rds.html">Clustered RDS &uarr;</a></li>
  </ol>
  <div class="pipeline-junction"><strong>Annotated cells</strong><a class="pipeline-resume" href="inputs/rds.html">Annotated RDS &larr;</a></div>
  <div class="pipeline-results">
    <div><a href="downstream/differential.html">Differential expression</a><p>Sample-level expression contrasts</p><div class="pipeline-followup"><span aria-hidden="true">&darr;</span><a href="downstream/enrichment.html">Functional enrichment</a><p>GO and KEGG over-representation</p></div></div>
    <div><a href="downstream/monocle3.html">Trajectory inference</a><p>Cell ordering and pseudotime</p></div>
    <div><a href="downstream/cellchat.html">Cell-cell communication</a><p>Ligand-receptor networks</p></div>
  </div>
  </div></section>
</div>

## Choose your starting point

<div class="starting-points" markdown="1">

[**FASTQ**<br>Paired-end reads requiring alignment and count generation.](inputs/fastq.md)

[**10x Matrix**<br>Existing count matrices ready for cell-level QC.](inputs/matrix.md)

[**Seurat RDS**<br>Resume after filtering, normalization, clustering, or annotation.](inputs/rds.md)

</div>

New users can begin with [Installation](installation.md) and [Quick Start](quickstart.md).
For an existing run, use [Configuration](reference/configuration.md),
[Workflow Targets](reference/targets.md), or [Troubleshooting](reference/troubleshooting.md).

## From broad labels to cell subtypes

SingleR reference labels are summarized at cluster level, then matched to eligible marker lineages.
UCell scoring evaluates candidate subtypes within each lineage. Supported subtype assignments
become final fine labels; insufficient evidence retains a broad label or `Unknown`.

The [Cell Annotation](workflow/annotation.md) page connects this procedure to the broad and fine
UMAP views, cluster-level score table, marker-expression heatmap, and composition reports.
These outputs provide complementary evidence for reviewing labels; finer labels alone do not
establish greater accuracy.

[![Original clusters, major cell types, and refined subtypes from the supplied annotation report.](assets/examples/annotation.png)](workflow/annotation.md#inspect-the-results)

**Annotation example.** The supplied results folder `5` shows the same embedding
colored by clusters, broad labels, and refined subtypes. Some labels remain broad or
`Unknown`; review the [annotation evidence](workflow/annotation.md#inspect-the-results)
before using the labels downstream.

## Inspect results before continuing

- **Before normalization:** Review cell and gene retention in the [filtering report](workflow/qc_filtering.md#inspect-the-results).
- **Before annotation:** Check batch representation and cluster structure in the [clustering preview](workflow/clustering.md#inspect-the-results).
- **Before biological interpretation:** Compare annotation evidence with the requirements of the chosen downstream analysis.

Stage pages explain what to inspect and which settings affect the result. The
[Output Reference](reference/outputs.md) provides file locations; example parameter values are
starting settings, not a validated choice for every dataset.
