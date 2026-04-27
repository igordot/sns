---
title: rna-snv
parent: Routes
---

# Route: rna-snv

Variant detection in RNA-seq data.
Can be run following `rna-star`.

Segments:

* Trim adapters and low quality bases (Trimmomatic).
* Align to the reference genome (STAR).
* Align to other species and common contaminants (fastq_screen).
* Generate normalized genome browser tracks.
* Remove duplicate reads (Sambamba).
* Realign and recalibrate (GATK).
* Determine depth of coverage (GATK).
* Call point mutations and small insertions/deletions (GATK HaplotypeCaller and LoFreq).

For somatic variant detection, follow with [wes-pairs-snv](wes-pairs-snv).

## Usage

Set up a new analysis (common across all routes).
If running for the first time, check the [detailed usage instructions](../usage) for an explanation of every step.

```
cd <project dir>
git clone --depth 1 https://github.com/igordot/sns
sns/generate-settings <genome>
sns/gather-fastqs <fastq dir>
```

Run `rna-snv` route.

```
sns/run rna-snv
```

Check for potential problems.

```
grep "ERROR:" logs-sbatch/*
```

## Output
