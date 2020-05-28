---
title: Usage
nav_order: 2
---

# Usage

## Brief summary

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

Check if the jobs are submitted and running.

```
squeue -u $USER
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

Specify the reference genome (`hg19` or `hg38` for human, `mm10` for mouse, `dm3` or `dm6` for fly).

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

This can be done while the pipeline is still running and should be done after the pipeline completes.
There should be no output if everything ran without problems.
If any errors are detected, examine the log file where they are found to see the full context.

If there is a problem with any of the results, delete the broken files and re-run SNS.
It will generate any missing output.
Similarly, you can add additional entries to the sample sheet and only the new ones will be processed when the route is re-run.

## Output

* Directories for different output types (such as BAMs or bigWigs) containing files for each sample.
* `summary-combined.*.csv`: Combined segment summaries table that provides a comprehensive overview of the project.
* `logs-*` directories: Most stdout/stderr output will be placed here. The information can be used for tracking progress and troubleshooting.

Each route has a description with more specific details.
