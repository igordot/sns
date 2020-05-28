---
title: RNA-seq demo
parent: Tutorials
nav_order: 1
---

# Seq-N-Slide RNA-seq demo

## About

This RNA-seq demo is based on data from tumor associated macrophages (TAMs) isolated from wildtype (WT) or Myeloid-specific Tet2 knockout (KO) mice ([GSE98964](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98964)).
The raw FASTQs have already been downloaded on BigPurple. Only a subset of reads (up to 10M) was kept to speed up the analysis.

## Prepare the pipeline

Create the project directory (choose any name):

```
mkdir proj_dir
```

Navigate to the created project directory:

```
cd proj_dir
```

Load `git` module (`git` is not available by default on BigPurple):

```
module add git
```

Download the pipeline code:

```
git clone --depth 1 https://github.com/igordot/sns
```

This will create an `sns` sub-directory in the current directory.

## Process the individual RNA-seq samples

Specify the reference genome (`mm10` for mouse):

```
sns/generate-settings mm10
```

This will create `settings.txt`, which contains the information about the different reference files.

Specify the location of the raw FASTQ files:

```
sns/gather-fastqs /gpfs/data/igorlab/tutorials/FASTQ-RNA
```

This will create `samples.fastq-raw.csv`, which contains the sample names and the corresponding files.
If a single sample (same sample name) has multiple FASTQs, each FASTQ (or FASTQ pair) will be listed on a separate line.
The FASTQs will be automatically merged.
Sample names are automatically detected based on the file names, but they can be edited in this file to be more readable.
All downstream file names will contain the sample names specified in this file.

On BigPurple, the raw sequencing data from GTC is usually deposited in `/gpfs/data/sequence/results/[lab]/[date]`.

Execute the pipeline (`rna-star` route for standard RNA-seq analysis):

```
sns/run rna-star
```

## Check progress

Check if the jobs are submitted and running:

```
squeue -u $USER
```

Check for errors:

```
grep -i "error:" logs-sbatch/*
```

The command will search all the log files for any errors. This can be done while the pipeline is still running and should be done after the pipeline completes. There should be no output if everything ran without problems. If any errors are detected, open the log file where they are found to see the full context.

## Perform differential expression analysis

After the pipeline is finished, edit `samples.groups.csv` and specify groups. The differential expression step will compare all groups against each other.

Run the differential expression step (`rna-star-groups-dge` route):

```
sns/run rna-star-groups-dge
```
