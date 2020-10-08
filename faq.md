---
title: FAQs
nav_order: 6
---

# FAQs

**Is it possible to use a custom reference genome?**

Yes.
First, try using a known genome like hg38 or mm10 to familiarize yourself with the workflow and the expected output.
There should not be any errors as long as a small fraction of reads successfully align.
Check `settings.txt` to see which reference files were actually used.
The pipeline uses the files that are listed there.
Each genome directory contains a README file with information about how the files were generated.
General notes and tips are available in the [reference genomes](https://github.com/igordot/reference-genomes) repository.
Run the pipeline again from a new project directory with a modified `settings.txt` that includes the relevant files.

**Why not use a "proper" workflow system (such as Snakemake, Nextflow, CWL, etc.)?**

Various workflow systems tend to be designed by developers for developers.
SNS is designed for users and is optimized for ease of use (as long as you are familiar with the linux/unix command-line interface):

* Simple to run: no installation, no setup, no dependencies.
This applies not only to the pipeline code, but also to the software tools and reference files.
* Simple to interpret: the output is organized and summarized to be human-readable.
The pipeline also includes custom checks to detect problems with the input data even if the tools ran successfully and generated valid output.
* Simple to troubleshoot: the executed commands are not hidden by wrappers and can be run independently of the pipeline.
* Simple to read and modify: written in standard bash, so learning new languages is not required.
