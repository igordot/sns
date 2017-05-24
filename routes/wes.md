# Route: wes

Alignment and variant detection for whole genome/exome/targeted sequencing data using BWA-MEM and GATK.

Segments:

* Trim adapters and low quality bases (Trimmomatic).
* Align to the reference genome (BWA-MEM).
* Remove duplicate reads (Sambamba).
* Realign and recalibrate (GATK).
* Determine fragment size distribution.
* Determine capture efficiency and depth of coverage (GATK).
* Call point mutations and small indels (GATK HaplotypeCaller and LoFreq).

## Usage

Navigate to a clean new project directory.

```
cd <project dir>
```

Download the code from GitHub.

```
git clone --depth 1 https://github.com/igordot/sns
```

Add a BED file of target regions to the project directory.

Generate a sample sheet of FASTQ files (`samples.fastq-raw.csv`).

```
sns/gather-fastqs <fastq dir>
```

Specify a reference genome, such as `hg19` or `mm10` (stored in `settings.txt`).

```
sns/generate-settings <genome>
```

Run `wes` route.

```
sns/run wes
```

Check for potential problems.

```
grep "ERROR:" logs-qsub/*
```

## Output

Results:

* `BAM-GATK-RA-RC`: Final BAM files (deduplicated, realigned, and recalibrated). Can be used for visual inspection of variants or additional analysis.
* `VCF-*`: VCF files for GATK HaplotypeCaller and LoFreq variant callers.
* `VCF-*-annot.all.txt`: Table of functionally annotated variants.
* `VCF-*-annot.coding.txt`: Table of coding region variants.
* `VCF-*-annot.nonsyn.txt`: Table of non-synonymous, frameshift, and splicing variants.
* `samples.pairs.csv`: Sample sheet for additional paired analysis (if necessary).

Run metrics:

* `summary-combined.wes.csv`: Summary table that includes the number of reads, alignment rate, fraction of PCR duplicates, capture efficiency (enrichment in targeted regions), and depth/evenness of coverage.
* `summary.qc-fragment-sizes.png`: Distribution of fragment sizes.
* `summary.VCF-*-annot.csv`: Total number of mutations for different variant callers.

Additional output (can usually be deleted or used for troubleshooting):

* `logs-*`: Logs and intermediate files for various segments.
* `samples.*.csv`: Sample sheet for segments that generate large files. The route will not attempt to generate the files listed. If the files were deleted to save space, additional samples can be added to the same analysis without reprocessing the older samples.
* `summary`: Summary files for individual samples and segments.
* `summary.*.csv`: Combined summary files for each segment.
* `QC-*`: Results of QC steps for individual samples.
* `FASTQ-CLEAN`: Merged FASTQs (one per sample).
* `FASTQ-TRIMMED`: Quality and adapter trimmed FASTQs.
* `BAM-*`: Intermediate BAM files.

General pipeline info: https://github.com/igordot/sns
