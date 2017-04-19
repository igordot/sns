# Seq-N-Slide: illumina sequencing data analysis pipelines

## Usage overview

Navigate to a clean new project directory. This is where all the results will end up.
```
cd <project dir>
```

Download the code from GitHub, which will create the `sns` sub-directory with all the code.
```
git clone --depth 1 https://github.com/igordot/sns
```

Scan a directory that contains FASTQ files to be used as input.
This can be run multiple times if there are FASTQs in different directories.
```
sns/gather-fastqs <fastq dir>
```
All found files will be added to the `samples.fastq-raw.csv` file, which can be modified to adjust sample names or remove samples.
The first column is the sample name.
The second column is the R1 FASTQ.
The third column is the R2 FASTQ (if available).
Each line contains a single FASTQ (or FASTQ pair for paired-end experiments).
If one sample has multiple FASTQs, each one will be on a different line.
Multiple FASTQs for the same sample will be merged based on sample name.

Specify a genome (only `hg19/mm10/dm3/dm6` are currently guaranteed to work).
```
sns/generate-settings <genome>
```

Run the analysis using a specific route (a set of analysis steps).
```
sns/run <route>
```

Check for potential problems.
```
grep "ERROR:" logs-qsub/*
```
There should be no matches if everything is okay.
If there are any results, check the specific log files where the errors are found for more info.

## Routes

Routes are different analysis workflows. Generic routes are sample-centric (same analysis is performed for each sample). Available routes:
* `rna-star`: RNA-seq using STAR. Generates BAMs, normalized bigWigs, counts matrix, and various QC metrics.
* `rna-rsem`: RNA-seq using RSEM. Generates FPKM/TPM/counts matrix and various QC metrics.
* `rna-snv`: RNA-seq variant detection using STAR and GATK. Generates BAMs, VCFs, and various QC metrics.
* `wgbs`: WGBS using Bismark.
* `rrbs`: RRBS using Bismark.
* `wes`: Whole genome/exome/targeted variant detection using BWA-MEM and GATK. Generates BAMs, VCFs, and various QC metrics.
* `atac`: ATAC-seq using Bowtie and MACS. Generates BAMs, bigWigs, peaks, nucleosome positions, and various QC metrics.
* `species`: Species/metagenomics/contamination analysis using Centrifuge with NCBI BLAST nt/nr nucleotide collection.

There are additional routes for comparing groups of samples after individual samples are processed with a generic route.
They depend on the output of the generic routes and must be run from the same directory.
Before running, manually add proper group names or pairs to the `samples.groups.csv` or `samples.pairs.csv` files (depending on the comparison type).
Available comparison routes:
* `rna-star-groups-dge`: Differential gene expression using DESeq2 for the `rna-star` results.
* `wes-pairs-snv`: Somatic variant detection for the `wes` results.

## Output

* Directories for different output types (such as BAMs or bigWigs) containing files for each sample.
* `summary.*.csv`: Most segments will generate a results summary file. This lets you know if the step completed and some relevant statistics.
* `summary-combined.*.csv`: Combined segment summaries. This table provides a comprehensive overview of the project and should reveal any potential problems.
* `samples.*.csv`: Most segments that generate large files will also generate a separate segment-specific sample sheet. If the files referenced in the sample sheet are missing, the route will not attempt to generate them. This can be useful if the files were deleted to save space, but you want to add more samples to the same analysis without reprocessing the older samples.
* `logs-*` directories: Most stdout/stderr output will be placed here. The information can be used for tracking progress and troubleshooting, but is generally useless.

## About

SNS is designed to work on NYULMC HPC cluster using the Sun Grid Engine job scheduler.
It may require significant modifications to work in other environments.

SNS consists of multiple routes (or workflows).
Each route contains multiple segments (or steps).

If there is a problem with any of the results, delete the broken files and re-run SNS.
It will generate any missing output.

Most output and sample sheets are in a CSV format for macOS Quick Look (spacebar file preview) compatibility.

There is a copy of the code in each project directory for reproducibility.
If you modify the code, the changes will not affect other projects.
If you repeat the analysis with more samples in the future, same code will be used.

## FAQs

Coming soon.

