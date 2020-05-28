# Seq-N-Slide: sequencing data analysis pipelines

Automated workflows for common sequencing-based (Illumina) protocols, such as RNA-seq, ChIP-seq, ATAC-seq, WGBS/RRBS methylation, whole genome/exome/targeted variant detection, and contaminant screening.

## Brief usage overview

Download the code.

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

For more details, check the full documentation at: https://igordot.github.io/sns
