# chum

[![Install with bioconda](https://img.shields.io/badge/Install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/chum/README.html)
[![Anaconda Version](https://anaconda.org/bioconda/chum/badges/version.svg)](http://bioconda.github.io/recipes/chum/README.html)
[![Build Status](https://github.com/fg-labs/chum/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/fg-labs/chum/actions/workflows/rust.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/fg-labs/chum/badge.svg?branch=main)](https://coveralls.io/github/fg-labs/chum?branch=main)
[![Language](https://img.shields.io/badge/language-rust-dea588.svg)](https://www.rust-lang.org/)

Evaluate the effectiveness of baits in a hybrid selection panel.

![Kauaʻi](.github/img/cover.jpg)

Install with mamba, conda, or run directly with pixi:

```bash
pixi exec \
    -c conda-forge -c bioconda \
    chum --help
```

## Introduction

Hybrid selection (also known as capture sequencing) relies on short oligonucleotide baits to enrich target regions from a sequencing library.
The efficiency of a hybrid selection panel depends on many bait-level properties such as nucleotide composition, sequence complexity, genomic uniqueness, overlap with repetitive elements, and thermodynamic stability of RNA/DNA secondary structures at hybridization temperature.
This toolkit provides subcommands for common tasks when managing hybrid selection experiments including the evaluation of bait performance.

## Scoring Baits (_In Silico_)

The subcommand `score` computes a comprehensive set of _in silico_ per-bait quality metrics and aggregates them per-target to help users evaluate, filter, and compare baits before committing to synthesis.
The primary input is a set of baits in BED, Interval List, or FASTA format.
BED and Interval List inputs require a reference FASTA to extract bait sequences for sequence-based metrics.
FASTA bait input works best if you embed the position-formatted genomic coordinates for the bait (`chrom:start-end`) as a comment in the FASTA record header.
The primary output is a tab-separated metrics file with one row per bait in the same order as the input.

#### Example Usage

Evaluate baits from a FASTA file against a reference genome with all optional analyses enabled:

```bash
❯ chum score \
    --baits "baits.fa" \
    --targets "targets.bed" \
    --reference "hs38GIABv3.fa" \
    --blast-db hs38GIABv3 \
    --rep-base "hs38GIABv3_rmsk.bed.gz" \
    --mappability "k36.umap.bedgraph.gz" \
    --oligo-fold \
    --per-bait "per-bait.tsv" \
    --per-target "per-target.tsv"
```

#### Features

- Accepts baits in BED, Interval List, or FASTA format (coordinates can be embedded in FASTA headers)
- Computes per-bait QC metrics: sequence content, secondary structure, mappability, BLAST specificity
- Multi-factor `bait_score` from 0.0 to 1.0 which may prompt discarding or redesigning of baits
- Per-target group aggregation (min/mean/max per numeric field) with `--targets` & `--per-target`
- Parallelizes bait evaluation across worker threads with `--threads` for very fast runtimes

###### Target Evaluation

When a target interval file is supplied (`--targets`), `chum` outputs per-target aggregate metrics including minimum, mean, and maximum of each numeric bait field across all overlapping (or nearby) baits.
Targets can be optionally padded which helps with linking near-baits to targets on either side of the original target interval with `--target-padding`.
The aggregate output file is specified with `--per-target`, which must be supplied whenever `--targets` is set.

When interpreting mean metric values, note that baits with no value for a field are excluded from the mean, which may skew results.
For example, if a target has three baits with `blast_hits` of `None`, `Some(2)`, and `None`, the mean `blast_hits` is reported as `2`, not `0.67`.

###### BLAST Alignment

When `blastn` is on the system path and a [BLAST](https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/) database is supplied (`--blast-db`), each bait sequence is aligned to the database and summary statistics from the top two BLAST hits are included in the output.
BLAST metrics help assess how specific a bait is to its intended target in the genome.

A high BLAST hit count does not necessarily indicate poor specificity.
A user must examine the percent identity and e-value of the first and second hits to judge whether off-target baiting is a concern.
Baits from FASTA input without embedded coordinates will have mappability and RepBase metrics computed from their BLAST top-hit coordinates when available.

###### BLAST Query Complexity Masking

BLAST uses the [DUST program](https://meme-suite.org/meme/doc/dust.html) by default to mask low-complexity regions in query sequences before seeding alignments, which limits statistically valid but biologically uninteresting hits.
Because `chum` evaluates baits against the whole genome including repeat-rich and low-complexity loci, DUST masking is disabled by default.
Enable it with `--blast-dust` when you want to suppress alignments driven purely by sequence composition.

###### Mappability

When a block-compressed and tabix-indexed bedGraph of mappability scores is supplied (`--mappability`), per-base mappability statistics are included in the output.
Mappability reflects how uniquely short sub-sequences of the bait map back to the reference genome and predicts how well sequenced fragments from the target region can be unambiguously placed.
It is recommended to use [Umap](https://bismap.hoffmanlab.org/) bedGraph files derived from multi-read measure and a _k_-mer size of 36.

###### Repeat Masker / RepBase

When a block-compressed and tabix-indexed BED file of [RepBase features](https://repeatbrowser.ucsc.edu/data/) is supplied (`--rep-base`), each bait is intersected with the feature set and overlapping repeat names are included in the output.
Repeat overlaps are important for predicting bait synthesis efficacy and the likelihood of alignment artifacts.
For example, simple repeats such as `(GGGGGC)n` can reduce hybrid selection efficiency, while elements like Alu or LINE/L1 can confound captures due to high endogenous copy number.

###### Oligo Secondary Structure

When `--oligo-fold` is set and `RNAfold` ([ViennaRNA](https://www.tbi.univie.ac.at/RNA/)) is on the system path, each bait sequence is folded at the temperature specified by `--oligo-fold-temp` (default 65 °C).
The minimum free energy (MFE, kcal/mol) and dot-bracket secondary structure string are included in the output.
A highly stable secondary structure (very negative MFE) may reduce hybridization efficiency.
Although `RNAFold` was originally designed to fold longer RNA oligonucleotides, `chum` is parameterized to use a DNA configuration by default for ssDNA baits.

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
