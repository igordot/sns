---
title: Output
nav_order: 3
---

# Output

Most segments will generate a sub-directory of results containing an output file for each sample.

Some segments may generate a logs sub-directory with stdout/stderr output or intermediate files which can be used for tracking progress and troubleshooting.
If the analysis is completed successfully, the logs directories can usually be deleted.

Standard output for all routes:

* `summary-combined.[route].csv`: Combined segment summary table that provides a comprehensive overview of the project.
* `summary.[segment].csv`: Summary table with important metrics for specific segments.
* `samples.[segment].csv`: Sample sheet for segments that generate large files. The route will not attempt to generate the files listed. If the files were deleted to save space, additional samples can be added to the same analysis without reprocessing the older samples.
* `FASTQ-CLEAN`: FASTQs that have been merged (if more than one per sample) and renamed to match the sample names in the sample sheet.

Common output:

* `FASTQ-TRIMMED`: Quality and adapter trimmed FASTQs.
* `BAM-*`: Original and post-processed BAM files.
* `QC-*`: Various QC metrics.

Each route has a description with details regarding its own output.
