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

For differential expression analysis, follow by running `rna-star-groups-dge`.

## Usage

Navigate to a clean new project directory.

```
cd <project dir>
```

Download the code from GitHub.

```
git clone --depth 1 https://github.com/igordot/sns
```

Generate a sample sheet of FASTQ files (`samples.fastq-raw.csv`).

```
sns/gather-fastqs <fastq dir>
```

Specify a reference genome, such as `hg19` or `mm10` (stored in `settings.txt`).

```
sns/generate-settings <genome>
```

Run `rna-star` route.

```
sns/run rna-star
```

Check for potential problems.

```
grep "ERROR:" logs-qsub/*
```

## Output

Results:

* `BAM-STAR`: BAM files. Can be used for visual inspection of individual reads or additional analysis.
* `BIGWIG`: BigWig files normalized to the total number of reads. Can be used for visual inspection of relative expression levels.
* `quant.featurecounts.counts.txt`: Matrix of raw counts for all genes and samples.

Run metrics:

* `summary-combined.rna-star.csv`: Summary table that includes the number of reads, unique and multi-mapping alignment rate, number of counts assigned to genes, fraction of coding/UTR/intronic/intergenic bases.
* `summary.fastqscreen.png`: Alignment rates for common species and contaminants.
* `summary.qc-picard-rnaseqmetrics.png`: Distribution of the bases within the transcripts to determine potential 5'/3' biases.

Additional output (can usually be deleted or used for troubleshooting):

* `logs-*`: Logs and intermediate files for various segments.
* `samples.*.csv`: Sample sheet for segments that generate large files. The route will not attempt to generate the files listed. If the files were deleted to save space, additional samples can be added to the same analysis without reprocessing the older samples.
* `summary`: Summary files for individual samples and segments.
* `summary.*.csv`: Combined summary files for each segment.
* `QC-*`: Results of QC steps for individual samples.
* `FASTQ-CLEAN`: Merged FASTQs (one per sample).
* `genes.featurecounts.txt`: Table of genes based on the reference GTF.
* `quant-*`: Raw counts for all genes for individual samples.
