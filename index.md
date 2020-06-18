---
title: Home
nav_order: 1
---

# Seq-N-Slide: sequencing data analysis pipeline

Seq-N-Slide is a set of automated workflows for common sequencing-based (Illumina) protocols, such as RNA-seq, ChIP-seq, ATAC-seq, WGBS/RRBS methylation, whole genome/exome/targeted variant detection, and contaminant screening.

SNS consists of multiple routes (or workflows).
Each route contains multiple segments (or steps).

Most output and sample sheets are in a CSV format for macOS Quick Look (spacebar file preview) compatibility.

Each project directory can have its own copy of the code in for reproducibility.
If you modify the code for one project, the changes will not affect other projects.
If you repeat the analysis with additional samples at a different time, same code will be used.

SNS is designed to work on the NYULMC HPC BigPurple cluster using the SLURM job scheduler.
It was originally made for the Phoenix cluster with the Sun Grid Engine job scheduler and that version has been archived on the [`phoenix` branch](https://github.com/igordot/sns/tree/phoenix).
It can work in other environments with modifications.

Source code is available at [github.com/igordot/sns](https://github.com/igordot/sns).
