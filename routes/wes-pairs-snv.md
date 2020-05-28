---
title: wes-pairs-snv
parent: Routes
---

# Route: wes-pairs-snv

Somatic point mutations and small indels for the `wes` results.

## Usage

After individual samples are processed with the `wes` route,
manually define tumor-normal pairs in the `samples.pairs.csv` sample sheet.

Run `wes-pairs-snv` route from the same directory as `wes`.

```
sns/run wes-pairs-snv
```

## Output

Results:

* `VCF-MuTect2`: VCF files for MuTect2 variant caller.
* `VCF-MuTect2-annot.all.txt`: Table of functionally annotated variants.
* `VCF-MuTect2-annot.coding.txt`: Table of coding region variants.
* `VCF-MuTect2-annot.nonsyn.txt`: Table of non-synonymous, frameshift, and splicing variants.

Run metrics:

* `summary.VCF-MuTect2-annot.csv`: Total number of mutations.
