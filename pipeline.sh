#!/bin/bash

# ==============================================================================
# FASTQ to VCF pipeline - Escherichia coli K-12 MG1655
#
# A variant calling pipeline: from raw sequencing reads (FASTQ) to variants
# (VCF), aligning a sample against a reference genome.
#
# Sample: E. coli, 50,000 paired-end reads (150 bp), subset to run SRR957824.
# Reference: E. coli K-12 MG1655 (NC_000913.3, NCBI RefSeq).
#
# Environment: WSL2 Ubuntu, conda environment 'fastq2vcf'
# Tools: fastqc, fastp (installed via bioconda)
#
# Author: Emmanuel Humberto
# ==============================================================================

# Note: data files (FASTQ, reference) are NOT versioned (see .gitignore).
# They can be re-downloaded following the steps in the README.

# ------------------------------------------------------------------------------
# Step 1: Quality control of raw reads (FastQC)
# ------------------------------------------------------------------------------
# Inspect the quality of the raw paired-end reads before any processing.
fastqc data/ecoli_R1.fastq.gz data/ecoli_R2.fastq.gz -o qc/

# Diagnosis: good quality in the body of the reads, with two expected technical
# residues at the ends -> quality dropping at the 3' end (normal Illumina decay)
# and TruSeq adapters (fragments shorter than read length). Both are treatable
# by trimming, not signs of bad data.

# ------------------------------------------------------------------------------
# Step 2: Read trimming (fastp)
# ------------------------------------------------------------------------------
# Remove adapters and trim low-quality 3' ends. Paired-end aware: R1 and R2 are
# processed together so mates stay synchronised.
fastp \
  -i data/ecoli_R1.fastq.gz \
  -I data/ecoli_R2.fastq.gz \
  -o data/ecoli_R1.trimmed.fastq.gz \
  -O data/ecoli_R2.trimmed.fastq.gz \
  --detect_adapter_for_pe \
  --cut_tail \
  --cut_tail_mean_quality 20 \
  --length_required 50 \
  --html qc/fastp_report.html \
  --json qc/fastp_report.json

# Result: 50,000 -> 42,148 read pairs kept (~15% removed, dominated by adapter
# dimers). Q30 rose from 86.5% to 95.3% (R1) and 75.0% to 92.3% (R2).

# ------------------------------------------------------------------------------
# Step 3: Quality control of trimmed reads (FastQC) - verification
# ------------------------------------------------------------------------------
# Confirm the trimming worked: the three red flags (per-base quality, adapter
# content, overrepresented sequences) all turn green. The only new warning is
# variable read length, which is the expected signature of trimming.
fastqc data/ecoli_R1.trimmed.fastq.gz data/ecoli_R2.trimmed.fastq.gz -o qc/

# ------------------------------------------------------------------------------
# Next steps (to be added): alignment to reference, variant calling -> VCF
# ------------------------------------------------------------------------------
