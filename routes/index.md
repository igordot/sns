---
title: Routes
nav_order: 4
has_children: true
has_toc: false
---

# Routes

Routes are the different analysis workflows.

## Basic routes

Generic routes are sample-centric where the same analysis is performed independently for each sample.
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

Some protocols require analysis of multiple samples simultaneously, such as tumor-normal or treatment-control combinations.
There are additional routes for comparing groups of samples after individual samples are processed with a generic sample-centric route.
These routes depend on the output of the generic routes and must be run from the same directory.
Available comparison routes:

* [rna-star-groups-dge](rna-star-groups-dge): Differential gene expression for the `rna-star` results.
* [wes-pairs-snv](wes-pairs-snv): Somatic variant detection for the `wes` or `rna-snv` results.
* [chip-pairs-peaks](chip-pairs-peaks): Peak calling for the `chip` results.

Each route has a description with more details regarding its output.
