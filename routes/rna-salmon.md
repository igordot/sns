---
title: rna-salmon
parent: Routes
---

# Route: rna-salmon

Transcript quantification of RNA-seq data using Salmon.

Segments:

* Trim adapters and low quality bases (Trimmomatic).
* Quantify (Salmon).
* Summarize estimates for and gene-level analysis (tximport).

## Usage

Set up a new analysis (common across all routes).

```
cd <project dir>
git clone --depth 1 https://github.com/igordot/sns
sns/generate-settings <genome>
sns/gather-fastqs <fastq dir>
```

Run `rna-salmon` route.

```
sns/run rna-star
```

Check for potential problems.

```
grep "ERROR:" logs-qsub/*
```

## Output

Primary results:

* `quant.salmon.tximport.counts.txt`: Matrix of estimated counts for all genes and samples.
* `quant.salmon.tximport.counts.txt`: Matrix of TPMs for all genes and samples.

Run metrics:

* `summary-combined.rna-salmon.csv`: Summary table that includes the number of reads, the mapping rate, and the number of detected genes.
