---
title: Routes
nav_order: 4
has_children: true
has_toc: false
---

# Routes

Routes are the different analysis workflows.

## Basic routes

Generic routes are sample-centric (same analysis is performed for each sample).
The available routes:

* [rna-star](rna-star): RNA-seq using STAR. Generates BAMs, normalized bigWigs, counts matrix, and various QC metrics.
* [rna-salmon](rna-salmon): RNA-seq using Salmon. Generates TPM/counts matrix and various QC metrics.
* [rna-rsem](rna-rsem): RNA-seq using RSEM. Generates FPKM/TPM/counts matrix and various QC metrics.
* [rna-snv](rna-snv): RNA-seq variant detection. Generates BAMs, VCFs, and various QC metrics.
* [chip](chip): ChIP-seq. Generates BAMs, bigWigs, peaks, and various QC metrics.
* [atac](atac): ATAC-seq. Generates BAMs, bigWigs, peaks, nucleosome positions, and various QC metrics.
* [wes](wes): Whole genome/exome/targeted variant detection. Generates BAMs, VCFs, and various QC metrics.
* [wgbs](rrbs): WGBS methylation analysis.
* [rrbs](rrbs): RRBS methylation analysis.
* [species](species): Species/metagenomics/contamination analysis.

## Routes for comparing groups of samples

There are additional routes for comparing groups of samples after individual samples are processed with a generic route.
They depend on the output of the generic routes and must be run from the same directory.
Before running, manually add proper group names or pairs to the `samples.groups.csv` or `samples.pairs.csv` files (depending on the comparison type).
Available comparison routes:

* [rna-star-groups-dge](rna-star-groups-dge): Differential gene expression using DESeq2 for the `rna-star` results.
* [wes-pairs-snv](wes-pairs-snv): Somatic variant detection for the `wes` results.
* [chip-pairs-peaks](chip-pairs-peaks): Peak calling for the `chip` results.

Each route has a description with more specific details.
