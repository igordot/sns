---
title: wes-pairs-cnv
parent: Routes
---

# Route: wes-pairs-cnv

Somatic copy number alterations for the `wes` results.

## Usage

After individual samples are processed with the `wes` route, manually define tumor-normal sample pairs in the `samples.pairs.csv` sample sheet.

Add a BED file defining the genomic regions of the probes/baits to the project directory.
The baited regions (or capture targets) are the loci your kit actually captures and should include flanking regions on either side of the targeted regions.

Run `wes-pairs-cnv` route from the same directory as `wes`.

```
sns/run wes-pairs-cnv
```
