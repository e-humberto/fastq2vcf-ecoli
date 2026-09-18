# FASTQ to VCF — Variant Calling Pipeline (*E. coli*)

A complete variant-calling pipeline: from raw sequencing reads to a VCF file.
Starting from raw paired-end reads of an *Escherichia coli* sample, this project
runs the full workflow (quality control, trimming, alignment, and variant
calling) to identify where the sample differs from a reference genome.

## Objective

Given raw sequencing reads from a bacterial sample, produce a reliable VCF file
documenting its genomic variants relative to a reference, and document every
quality decision along the way.

## Data

- **Sample:** *E. coli* O157 (a pathogenic strain), 50,000 paired-end reads
  (150 bp), a subset of SRA run SRR957824.
- **Reference:** *E. coli* K-12 MG1655 (NCBI RefSeq `NC_000913.3`).

The sample and the reference are **different strains of the same species**. This
was a deliberate choice, and it shapes the results in an instructive way (see
Results).

## Pipeline

| Step | Tool | Purpose |
|------|------|---------|
| 1. Quality control | FastQC | Inspect raw read quality |
| 2. Trimming | fastp | Remove adapters and low-quality read ends |
| 3. QC verification | FastQC | Confirm trimming worked |
| 4. Reference indexing | BWA | Prepare the reference for alignment |
| 5. Alignment | BWA-MEM | Map reads to the reference (SAM) |
| 6. Sort & index | samtools | Coordinate-sorted, indexed BAM |
| 7. Variant calling | bcftools | Produce the VCF |

The full commands are in [`pipeline.sh`](pipeline.sh).

## Results

- **Trimming:** 50,000 → 42,148 read pairs retained (~15% removed, mostly
  adapter dimers). Q30 rose from 86.5% to 95.3% (R1) and 75.0% to 92.3% (R2).
- **Alignment:** 77.3% of reads mapped to the reference.
- **Variants:** 52,814 raw variants called.

### Reading the results honestly

The **77.3% mapping rate** and the **large number of variants** tell the same
story: the sample (O157) and the reference (K-12) are related but distinct
strains. Reads from O157-specific regions have no counterpart in the K-12
reference and therefore do not map; and even in shared regions, the two strains
differ at many positions. Both figures reflect real biological divergence, not a
pipeline error.

The 52,814 variants are the **raw, unfiltered** call set. A rigorous downstream
analysis would filter by quality (`QUAL`) and read depth (`DP`) to retain the
confident subset, many raw calls rest on very low coverage (e.g. `DP=1`) and
would not survive filtering.

## Tools

FastQC · fastp · BWA · samtools · bcftools

## Notes

This is a portfolio project built to demonstrate the end-to-end variant-calling
workflow with transparent, documented decisions at each step. Data files (FASTQ,
reference, BAM, VCF) are not versioned here; the pipeline reproduces them from
the sources described above.
