---
title: chip
parent: Routes
nav_order: 5
---

# Route: chip

ChIP-seq analysis.

Segments:

* Align to the reference genome (Bowtie2).
* Remove duplicate reads (Sambamba).
* Generate genome browser tracks.

## Usage

Navigate to a clean new project directory, download the code, specify a reference genome, and generate a sample sheet.
These steps are described in more detail on the [main page](https://github.com/igordot/sns).

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

Results:

* `BAM-DD`: Deduplicated BAM files. Can be used for additional analysis.
* `BIGWIG`: BigWig files normalized to the total number of reads. Can be used for visual inspection of peaks.
* `samples.pairs.csv`: Sample sheet for defining treatment-control pairs for peak calling.

Run metrics:

* `summary-combined.chip.csv`: Summary table that includes the number of reads, alignment rate, and fraction of PCR duplicates.

Additional output (can usually be deleted or used for troubleshooting):

* `logs-*`: Logs and intermediate files for various segments.
* `samples.*.csv`: Sample sheet for segments that generate large files. The route will not attempt to generate the files listed. If the files were deleted to save space, additional samples can be added to the same analysis without reprocessing the older samples.
* `summary`: Summary files for individual samples and segments.
* `summary.*.csv`: Combined summary files for each segment.
* `FASTQ-CLEAN`: Merged FASTQs (one per sample).
* `BAM`: Initial non-deduplicated BAM files.
