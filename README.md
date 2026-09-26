# Long-read RNA-seq Nextflow workflow

This repository contains a Nextflow workflow developed for DSA4262 Assignment 1. It processes Oxford Nanopore long-read RNA-sequencing data using Minimap2, Samtools and Bambu.

The workflow was adapted from Jonathan Göke's workshop workflow and extended to support multiple samples, protocol-specific alignment parameters, BAM sorting and indexing, alignment quality control, and Bambu runs with or without reference annotations.

## Workflow

For each sample, the workflow performs:

1. Splice-aware alignment to the reference genome using Minimap2.
2. SAM-to-BAM conversion, coordinate sorting and indexing using Samtools.
3. Alignment quality control using `samtools quickcheck`, `flagstat` and `stats`.
4. Joint transcript discovery and quantification using Bambu.

## Required inputs

- Oxford Nanopore FASTQ files
- Reference genome in FASTA format
- Reference transcript annotation in GTF format
- CSV samplesheet containing `sample_id`, `protocol` and `fastq`

The `protocol` column must contain either `directRNA` or `cDNA`. Direct-RNA reads are aligned using `-ax splice -uf -k14`, while cDNA reads are aligned using `-ax splice`.

An example is provided in `samplesheet.example.csv`.

## Software

- Nextflow 26.04.6
- Java 21
- Minimap2 2.26-r1175
- Samtools 1.19.2
- R 4.3.3
- Bambu 3.4.1
- Bioconductor 3.18

## Running with annotations

```bash
nextflow run main.nf \
    -with-report report_with_annotations.html \
    --samples /absolute/path/to/samplesheet.csv \
    --refFa /absolute/path/to/reference.fa \
    --refGtf /absolute/path/to/reference.gtf \
    --outdir /absolute/path/to/results \
    --bambu_mode with_annotations
```

## Running without annotations using cached processes

Run this command from the same directory after completing the annotated run:

```bash
nextflow run main.nf \
    -resume \
    -with-report report_without_annotations.html \
    --samples /absolute/path/to/samplesheet.csv \
    --refFa /absolute/path/to/reference.fa \
    --refGtf /absolute/path/to/reference.gtf \
    --outdir /absolute/path/to/results \
    --bambu_mode without_annotations
```

The `-resume` option allows Nextflow to reuse unchanged alignment, BAM conversion and QC processes from its cache. Only the modified Bambu process is executed again.

## Outputs

The workflow produces:

- Sorted and indexed BAM files
- `samtools quickcheck`, `flagstat` and `stats` QC files
- Extended transcript annotations in GTF format
- Transcript-level read counts
- Gene-level read counts
- A saved Bambu `RangedSummarizedExperiment` object

## Important note

The workflow currently uses `NDR = 1` for Bambu transcript discovery. This is permissive and was retained from the adapted workshop implementation. A lower or automatically determined NDR should be considered for analyses where stronger control of false-positive novel transcripts is required.
