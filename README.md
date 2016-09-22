# Seq-N-Slide: illumina sequencing data analysis pipelines

## Usage overview

Navigate to a clean new project directory.
```
cd <project dir>
```
This is where all the results will end up.

Download the code from the GitHub repository.
```
git clone --depth 1 https://github.com/igordot/sns
```
The project directory should now contain the `sns` sub-directory with all the code.

Find input FASTQ files in a given directory.
```
sns/gather-fastqs <fastq dir>
```
Can be run multiple times if there are multiple FASTQ directories. All found files will be added to the `samples.fastq-raw.csv` file, which can be modified to adjust sample names or remove samples. The first column is the sample name. The second column is the R1 FASTQ. The third column is the R2 FASTQ (for paired-end reads). Multiple FASTQs for the same sample will be merged.

Specify a genome (only `hg19/mm10/dm3/dm6` are currently guaranteed to work).
```
sns/generate-settings <genome>
```
This will create a `settings.txt` file.

Run the analysis using a specific route (a set of analysis steps).
```
sns/run <route>
```
Progress logs are available in the `logs-*` directories.

Check for any errors when it's done running or while it's still running.
```
grep -i "err" logs-qsub/*
```

## Routes

Routes are different analysis types. Generic routes are sample-based (same analysis is repeated for each sample). Available routes:
- `rna-star` - RNA-seq using STAR. Generates BAMs, normalized bigWigs, counts matrix, and various QC metrics.
- `rna-rsem` - RNA-seq using RSEM. Generates FPKM/TPM/counts matrix and various QC metrics.
- `rrbs` - RRBS using Bismark.
- `wgbs` - WGBS using Bismark.
- `wes` - Whole genome/exome/targeted variant discovery using BWA-MEM and GATK.
- `atac` - ATAC-seq using Bowtie and MACS. Generates BAMs, bigWigs, peaks, and various QC metrics.

There are also comparison routes for comparing multiple samples. Unlike sample names, group designations are difficult to determine automatically. You must modify the `samples.groups.csv` file with proper group names. Available comparison routes:
- `comp-rna-star-dge` - Differential gene expression using DESeq2 using the `rna-star` results.

## Output

- Directories for different output types (such as BAMs or bigWigs) containing files for each sample.
- `logs-*` directories - Most stdout/stderr output will be placed here. The information can be used for tracking progress and troubleshooting, but is generally useless.
- `summary.*.csv` - Most segments will generate a results summary file. This lets you know if the step completed and some relevant statistics.
- `summary-combined.*.csv` - Combined segment summaries.
- `samples.*.csv` - Most segments that generate large files will also generate a separate segment-specific sample sheet. If the files referenced in the sample sheet are missing, the route will not attempt to generate them. For example, this can be useful if the files were deleted to save space, but you want to add more samples to the same analysis without reprocessing the previous samples.

## About

SNS is designed to work on NYULMC HPC cluster using the Sun Grid Engine job scheduler. It may require significant modifications to work in other environments.

SNS consists of multiple routes (or workflows). Each route contains multiple segments (or steps).

If there is a problem with any of the results, delete the broken files and re-run SNS. It will generate any missing output.

Most output and sample sheets are in a CSV format for macOS Quick Look (spacebar file preview) compatibility.

There is a separate copy of the code in each project directory for reproducibility. If you modify the code, the changes will not affect other projects. If you repeat the analysis with more samples in the future, same code will be used.

## FAQs

Coming soon.

