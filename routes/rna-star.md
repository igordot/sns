---
title: rna-star
parent: Routes
---

# Route: rna-star

Alignment and quantification of RNA-seq data.

Segments:

* Trim adapters and low quality bases (Trimmomatic).
* Align to the reference genome (STAR).
* Align to other species and common contaminants (fastq_screen).
* Generate normalized genome browser tracks.
* Determine the distribution of the bases within the transcripts and 5'/3' biases (Picard).
* Determine if the library is stranded and the strand orientation.
* Generate genes-samples counts matrix (featureCounts).

For differential expression analysis, follow with [rna-star-groups-dge](rna-star-groups-dge).

## Usage

Set up a new analysis (common across all routes).
If running for the first time, check the [detailed usage instructions](../usage) for an explanation of every step.

```
cd <project dir>
git clone --depth 1 https://github.com/igordot/sns
sns/generate-settings <genome>
sns/gather-fastqs <fastq dir>
```

Run `rna-star` route.

```
sns/run rna-star
```

Check for potential problems.

```
grep "ERROR:" logs-sbatch/*
```

## Output

Primary results:

* `BAM-STAR`: BAM files. Can be used for visual inspection of individual reads or additional analysis.
* `BIGWIG`: BigWig files normalized to the total number of reads. Can be used for visual inspection of relative expression levels.
* `quant.featurecounts.counts.txt`: Matrix of raw counts for all genes and samples.

Run metrics:

* `summary-combined.rna-star.csv`: Summary table that includes the number of reads, unique and multi-mapping alignment rate, number of counts assigned to genes, fraction of coding/UTR/intronic/intergenic bases.
* `summary.fastqscreen.png`: Alignment rates for common species and contaminants.
* `summary.qc-picard-rnaseqmetrics.png`: Distribution of the bases within the transcripts to determine potential 5'/3' biases.

Additional output (can usually be deleted or used for troubleshooting):

* `genes.featurecounts.txt`: Table of genes based on the reference GTF.
* `quant-*`: Raw counts for all genes for individual samples.
