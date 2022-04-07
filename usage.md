---
title: Usage
nav_order: 2
---

# Usage

## Brief overview

Download the code in an empty project directory.

```
git clone --depth 1 https://github.com/igordot/sns
```

Generate a sample sheet based on a directory of FASTQ files.

```
sns/gather-fastqs <fastq_dir>
```

Specify the reference genome.

```
sns/generate-settings <genome>
```

Run the analysis using a specific route.

```
sns/run <route>
```

Check for problems.

```
grep "ERROR:" logs-sbatch/*
```

## Detailed description for new users

Add git to the environment (git is not available on BigPurple by default).

```
module add git
```

Navigate to a clean new project directory. This is where all the results will end up.

```
cd <project_dir>
```

Download the code from GitHub, which will create the `sns` sub-directory with all the pipeline code.

```
git clone --depth 1 https://github.com/igordot/sns
```

Search a directory of FASTQ files to be used as input and generate a sample sheet.

```
sns/gather-fastqs <fastq_dir>
```

All found files will be added to the `samples.fastq-raw.csv` file.
This command can be run multiple times if there are FASTQs in different directories.
The file can be edited to change the sample names, remove samples, or manually add samples.
All downstream file names will contain the sample names specified in this file.
The first column is the sample name, the second column is the R1 FASTQ, and the third column is the R2 FASTQ (if available).
Each line contains a single FASTQ (or FASTQ pair for paired-end experiments).
If a single sample has multiple FASTQs, each one will be on a different line.
Multiple FASTQs for the same sample will be merged based on the sample name.

Specify the reference genome (such as `hg38` or `hg19` for human, `mm10` for mouse, `dm6` or `dm3` for fly).

```
sns/generate-settings <genome>
```

This will create `settings.txt`, which contains the information about the reference files and certain project settings.

Run the analysis using a specific route (a set of analysis steps).

```
sns/run <route>
```

Check if the jobs are submitted and running.

```
squeue -u $USER
```

Check for potential problems.

```
grep "ERROR:" logs-sbatch/*
```

Checking for errors can be started as soon as the pipeline starts running.
It needs to be done after all the jobs complete.
There should be no output from this command if everything ran without problems.
If any errors are detected, examine the full log file where they are found to see the full context.

If there is a problem with any of the results, delete the problematic files and re-run SNS.
The pipeline will skip existing files and will generate any missing output.
Similarly, you can add additional entries to the sample sheet and only the new ones will be processed when the route is re-run.
