---
title: FAQs
nav_order: 6
---

# FAQs

**Is it possible to use a custom reference genome?**

Yes.
Before using a custom reference, try using a commonly used genome like hg38 or mm10 to familiarize yourself with the workflow and the expected output.
Check `settings.txt` to see which reference files were actually used.
You can run the pipeline from a new project directory with a pre-filled `settings.txt` that includes custom reference files.

The `settings.txt` file is created by the `generate-settings` step and additional entries are automatically added as needed while the pipeline runs for the first time.
The pipeline will look for relevant reference files in the specified genome directory.
Each genome directory contains a README file with information about how the reference files were generated.
General notes and tips are available in the [reference genomes](https://github.com/igordot/reference-genomes) repository.

**Can I run SNS anywhere?**

SNS is designed to work on the NYULMC HPC UltraViolet (formerly BigPurple) cluster using the SLURM job scheduler.
It was originally made for the Phoenix cluster with the Sun Grid Engine job scheduler and that version has been archived on the [`phoenix` branch](https://github.com/igordot/sns/tree/phoenix).
It can be modified to work in other environments, but the changes may not be trivial.
The computational steps expect that certain modules are available and link to various pre-generated reference files.

**Why not use a "proper" workflow system (such as Snakemake, Nextflow, CWL, etc.)?**

Workflow systems tend to be designed by developers for developers.
SNS is designed with novice users in mind and optimized for ease of use (assuming some familiarity with the linux/unix command-line interface):

* Simple to run: no installation, no setup, no dependencies.
This applies not only to the pipeline code, but also to the software tools and reference files.
* Simple to interpret: the output is organized and summarized to be human-readable.
The pipeline also includes custom checks to detect problems with the input data beyond simply checking if the tools ran successfully and generated valid output.
* Simple to troubleshoot: the executed commands are not hidden by wrappers and can be run independently of the pipeline.
* Simple to read and modify: written in standard bash, so learning new languages is not required.

**What if I have other questions?**

Any issues or concerns can be posted on [GitHub](https://github.com/igordot/sns/issues).
