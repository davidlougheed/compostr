# trlib

`trlib` is a tandem repeat motif decomposition and consensus sequence calculation library, primarily for use by
[STRkit](https://github.com/davidlougheed/strkit/) and `telokit`.


## Motivation

Tandem repeats (TRs) are sequences of DNA composed of a motif (or possibly a set of different motifs) repeated several
times in a row. are a common pattern in genomes. TRs can be characterized by their repeat copy number, motif
composition, and motif "purity", and can be disease-causing when the repeat expands.

One of the most well-known examples of this is [Huntington's disease](https://strchive.org/loci/hd_htt/), is caused by
a repeat expansion in a
[region](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hub_3671779_hs1&position=chr4:3073604-3073687) that normally looks
like this:

```
(CAG)6-26
```

If this repeat expands to, or is inherited at, 36+ copies of CAG in one of the two copies of the locus, that individual
will develop Huntington's disease.

As an example of how both length and motif composition can matter, consider
[Spinocerebellar ataxia type 37](https://strchive.org/loci/sca37_dab1/), where pathogenicity comes from an alternate
pentanucleotide motif stretch `(GAAAT)n`, rather than an expansion of the "canonical" motif `(AAAAT)n`.

**Therefore**, it is important to be able to decompose a TR sequence into the motifs that compose it.

Doing this is complicated by sequencing errors (string substitutions/insertions/deletions) and biology being messy as
usual, meaning we may get imperfect or erroneous motif copies.

**The goal here** is to build a library which can perform this motif decomposition accurately, quickly, and in a manner
which allows for non-canonical motifs or even non-repeat DNA inserted into the sequence.


## Copyright Notice

&copy; David Lougheed 2026.

### Notice

This library is licensed under the terms of the [Lesser GNU Public License 3.0](./LICENSE).
