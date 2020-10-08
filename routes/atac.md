---
title: atac
parent: Routes
---

# Route: atac

ATAC-seq analysis.

Segments:

* Align to the reference genome (Bowtie2).
* Remove duplicate reads (Sambamba).
* Generate genome browser tracks.
* Call peaks (MACS and HMMRATAC).

## Usage

Set up a new analysis (common across all routes).

```
cd <project dir>
git clone --depth 1 https://github.com/igordot/sns
sns/generate-settings <genome>
sns/gather-fastqs <fastq dir>
```

Run `atac` route.

```
sns/run atac
```

Check for potential problems.

```
grep "ERROR:" logs-sbatch/*
```

## Output

Primary results:

* `BAM-DD`: Deduplicated BAM files. Can be used for additional analysis.
* `BIGWIG`: BigWig files normalized to the total number of reads. Can be used for visual inspection of peaks.
* `peaks-*`: Peaks.

Run metrics:

* `summary-combined.atac.csv`: Summary table that includes the number of reads, alignment rate, and fraction of PCR duplicates.
* `summary.qc-fragment-sizes.png`: Distribution of fragment sizes.
