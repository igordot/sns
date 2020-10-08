---
title: chip-pairs-peaks
parent: Routes
---

# Route: chip-pairs-peaks

Peak calling for the `chip` results.

## Usage

After individual samples are processed with the `chip` route, manually define the type of peak in `settings.txt` and treatment-control pairs in the `samples.pairs.csv` sample sheet.

ChIP-seq can be used to detect transcription factor binding sites (narrow peaks) as well as enriched regions associated with most histone modification (broad peaks).
According to the ENCODE standards, broad marks include H3F3A, H3K27me3, H3K36me3, H3K4me1, H3K79me2, H3K79me3, H3K9me1, H3K9me2, and H4K20me1.
Narrow marks include H2AFZ, H3ac, H3K27ac, H3K4me2, H3K4me3, and H3K9ac.

The sample sheet requires a treatment and a control (genomic input or mock IP) sample. If no control is available, specify the same sample in both columns.

Run `chip-pairs-peaks` route from the same directory as `chip`.

```
sns/run chip-pairs-peaks
```

## Output

Primary results:

* `peaks-*`: Peaks.
