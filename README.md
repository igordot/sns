# Seq-N-Slide: illumina sequencing data analysis pipelines

Automated workflows for common sequencing-based protocols, such as RNA-seq, ATAC-seq, WGBS/RRBS methylation, whole genome/exome/targeted variant detection, and contaminant screening.

## Usage overview

Add git to the environment (git is not available on BigPurple by default).

```
module add git
```

Navigate to a clean new project directory. This is where all the results will end up.

```
cd <project dir>
```

Download the code from GitHub, which will create the `sns` sub-directory with all the relevant files.

```
git clone --depth 1 https://github.com/igordot/sns
```

Scan a directory that contains FASTQ files to be used as input.
This command can be run multiple times if there are FASTQs in different directories.

```
sns/gather-fastqs <fastq dir>
```

All found files will be added to the `samples.fastq-raw.csv` file.
It can be modified to change sample names, remove samples, or manually add samples.
The first column is the sample name, the second column is the R1 FASTQ, and the third column is the R2 FASTQ (if available).
Each line contains a single FASTQ (or FASTQ pair for paired-end experiments).
If one sample has multiple FASTQs, each one will be on a different line.
Multiple FASTQs for the same sample will be merged based on sample name.

Specify the reference genome, such as `hg19`, `hg38`, `mm10`, `dm3`, or `dm6`.

```
sns/generate-settings <genome>
```

Run the analysis using a specific route (a set of analysis steps).

```
sns/run <route>
```

Check for potential problems.

```
grep "ERROR:" logs-sbatch/*
```

There should be no matches.
If there are any results, there was a problem.
Check the specific log files where the errors are found for more info.

## Routes

Routes are the different analysis workflows.
Generic routes are sample-centric (same analysis is performed for each sample).
The available routes:

* [rna-star](https://github.com/igordot/sns/blob/master/routes/rna-star.md): RNA-seq using STAR. Generates BAMs, normalized bigWigs, counts matrix, and various QC metrics.
* [rna-rsem](https://github.com/igordot/sns/blob/master/routes/rna-rsem.md): RNA-seq using RSEM. Generates FPKM/TPM/counts matrix and various QC metrics.
* [rna-salmon](https://github.com/igordot/sns/blob/master/routes/rna-salmon.md): RNA-seq using Salmon. Generates TPM/counts matrix and various QC metrics.
* [rna-snv](https://github.com/igordot/sns/blob/master/routes/rna-snv.md): RNA-seq variant detection. Generates BAMs, VCFs, and various QC metrics.
* [atac](https://github.com/igordot/sns/blob/master/routes/atac.md): ATAC-seq. Generates BAMs, bigWigs, peaks, nucleosome positions, and various QC metrics.
* [wes](https://github.com/igordot/sns/blob/master/routes/wes.md): Whole genome/exome/targeted variant detection. Generates BAMs, VCFs, and various QC metrics.
* [wgbs](https://github.com/igordot/sns/blob/master/routes/rrbs.md): WGBS methylation analysis.
* [rrbs](https://github.com/igordot/sns/blob/master/routes/rrbs.md): RRBS methylation analysis.
* [species](https://github.com/igordot/sns/blob/master/routes/species.md): Species/metagenomics/contamination analysis.

There are additional routes for comparing groups of samples after individual samples are processed with a generic route.
They depend on the output of the generic routes and must be run from the same directory.
Before running, manually add proper group names or pairs to the `samples.groups.csv` or `samples.pairs.csv` files (depending on the comparison type).
Available comparison routes:

* [rna-star-groups-dge](https://github.com/igordot/sns/blob/master/routes/rna-star-groups-dge.md): Differential gene expression using DESeq2 for the `rna-star` results.
* [wes-pairs-snv](https://github.com/igordot/sns/blob/master/routes/wes-pairs-snv.md): Somatic variant detection for the `wes` results.

## Output

* Directories for different output types (such as BAMs or bigWigs) containing files for each sample.
* `summary-combined.*.csv`: Combined segment summaries table that provides a comprehensive overview of the project.
* `logs-*` directories: Most stdout/stderr output will be placed here. The information can be used for tracking progress and troubleshooting.

Each route has a description with more specific details.

## About

SNS is designed to work on the NYULMC HPC BigPurple cluster using the SLURM job scheduler.
It was originally made for the Phoenix cluster with the Sun Grid Engine job scheduler and that version has been archived on the [`phoenix` branch](https://github.com/igordot/sns/tree/phoenix).
It may require significant modifications to work in other environments.

SNS consists of multiple routes (or workflows).
Each route contains multiple segments (or steps).

If there is a problem with any of the results, delete the broken files and re-run SNS.
It will generate any missing output.
Similarly, you can add additional entries to the sample sheet and only the new ones will be processed when the route is re-run.

Most output and sample sheets are in a CSV format for macOS Quick Look (spacebar file preview) compatibility.

There is a copy of the code in each project directory for reproducibility.
If you modify the code for one project, the changes will not affect other projects.
If you repeat the analysis with more samples in the future, same code will be used.

## FAQs

Coming soon.

