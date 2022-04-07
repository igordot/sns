---
title: chip
parent: Routes
---

# Route: chip

ChIP-seq analysis.

Segments:

* Align to the reference genome (Bowtie2).
* Remove duplicate reads (Sambamba).
* Generate genome browser tracks.

For peak calling, follow with [chip-pairs-peaks](chip-pairs-peaks).

## Usage

Set up a new analysis (common across all routes).
If running for the first time, check the [detailed usage instructions](../usage) for an explanation of every step.

```
cd <project dir>
git clone --depth 1 https://github.com/igordot/sns
sns/generate-settings <genome>
sns/gather-fastqs <fastq dir>
```

Run `chip` route.

```
sns/run chip
```

Check for potential problems.

```
grep "ERROR:" logs-sbatch/*
```

## Output

Primary results:

* `BAM-DD`: Deduplicated BAM files. Can be used for additional analysis.
* `BIGWIG`: BigWig files normalized to the total number of reads. Can be used for visual inspection of peaks.

Run metrics:

* `summary-combined.chip.csv`: Summary table that includes the number of reads, alignment rate, and fraction of PCR duplicates.
