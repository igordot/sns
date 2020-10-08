---
title: rrbs
parent: Routes
---

# Route: rrbs/wgbs

Alignment and methylation calling for reduced representation or whole genome bisulfite sequencing data.

Segments:

* For WGBS, trim adapters and low quality bases (Trimmomatic).
* For RRBS, trim adapters, low quality bases, and the cytosine artificially introduced in the end-repair step during the library preparation (Trim Galore).
* Align to the bisulfite converted reference genome (Bismark and Bowtie2).
* For WGBS, remove duplicate reads (Bismark).
* Extract methylation calls (Bismark).
* Generate graphical report for alignment, deduplication, and methylation extraction (Bismark).
* Generate methylation calls genome browser tracks.

## Usage

Set up a new analysis (common across all routes).

```
cd <project dir>
git clone --depth 1 https://github.com/igordot/sns
sns/generate-settings <genome>
sns/gather-fastqs <fastq dir>
```

Run `rrbs` or `wgbs` route.

```
sns/run rrbs
```

Check for potential problems.

```
grep "ERROR:" logs-qsub/*
```

## Output

Primary results:

* `BAM-Bismark`: Alignment files.
* `BIGWIG-Bismark`: Methylation ratio bigWig files. Can be used for visual inspection as genome browser tracks.
* `meth-Bismark-*`: Methylation calls.

Run metrics:
 
* `summary-combined.rrbs.csv`: Summary table that includes the number of reads, alignment rate, fraction of PCR duplicates, number of covered Cs, and Cs methylated in different contexts.
* `Bismark-report`: Graphical report for Bismark steps (alignment, deduplication, and methylation extraction).
