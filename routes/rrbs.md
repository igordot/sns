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

Run `rrbs` or `wgbs` route.

```
sns/run rrbs
```

Check for potential problems.

```
grep "ERROR:" logs-qsub/*
```

## Output

Results:

* `BAM-Bismark`: Alignment files.
* `BIGWIG-Bismark`: Methylation ratio bigWig files. Can be used for visual inspection as genome browser tracks.
* `meth-Bismark-*`: Methylation calls.

Run metrics:
 
* `summary-combined.rrbs.csv`: Summary table that includes the number of reads, alignment rate, fraction of PCR duplicates, number of covered Cs, and Cs methylated in different contexts.
* `Bismark-report`: Graphical report for Bismark steps (alignment, deduplication, and methylation extraction).

Additional output (can usually be deleted or used for troubleshooting):

* `logs-*`: Logs and intermediate files for various segments.
* `samples.*.csv`: Sample sheet for segments that generate large files. The route will not attempt to generate the files listed. If the files were deleted to save space, additional samples can be added to the same analysis without reprocessing the older samples.
* `summary`: Summary files for individual samples and segments.
* `summary.*.csv`: Combined summary files for each segment.
* `FASTQ-CLEAN`: Merged FASTQs (one per sample).
* `FASTQ-TRIMMED`: Quality and adapter trimmed FASTQs.

General pipeline info: https://github.com/igordot/sns
