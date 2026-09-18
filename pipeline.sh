#!/bin/bash

# ==============================================================================
# FASTQ to VCF - Variant Calling Pipeline
# Organism: Escherichia coli
#
# From raw paired-end sequencing reads to a VCF file, aligning a sample against
# a reference genome and calling the variants that distinguish them.
#
# Sample:    E. coli O157 (pathogenic strain), 50,000 paired-end reads (150 bp),
#            subset of SRA run SRR957824.
# Reference: E. coli K-12 MG1655 (NCBI RefSeq NC_000913.3).
#
# Tools: FastQC, fastp, BWA, samtools, bcftools.
#
# Note: sample and reference are different strains of the same species. This is
# intentional; it illustrates how strain divergence affects mapping rate and
# variant counts (documented in the README).
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Quality control of the raw reads
# ------------------------------------------------------------------------------
# Inspect read quality before any processing.
fastqc data/ecoli_R1.fastq.gz data/ecoli_R2.fastq.gz -o qc/

# Diagnosis: good quality in the body of the reads, with two expected technical
# residues at the read ends -> quality decay at the 3' end (normal for Illumina)
# and TruSeq adapters (fragments shorter than the read length). Both are
# treatable by trimming, not signs of bad data.


# ------------------------------------------------------------------------------
# 2. Read trimming
# ------------------------------------------------------------------------------
# Remove adapters and trim low-quality 3' ends. Paired-end aware: R1 and R2 are
# processed together so mates stay synchronised.
fastp -i data/ecoli_R1.fastq.gz -I data/ecoli_R2.fastq.gz \
      -o data/ecoli_R1.trimmed.fastq.gz -O data/ecoli_R2.trimmed.fastq.gz \
      --detect_adapter_for_pe --cut_tail --cut_tail_mean_quality 20 \
      --length_required 50 \
      --html qc/fastp_report.html --json qc/fastp_report.json

# Result: 50,000 -> 42,148 read pairs kept (~15% removed, dominated by adapter
# dimers). Q30 rose from 86.5% to 95.3% (R1) and 75.0% to 92.3% (R2).


# ------------------------------------------------------------------------------
# 3. Quality control of the trimmed reads (verification)
# ------------------------------------------------------------------------------
# Confirm the trimming worked: the three FastQC red flags turn green. The only
# new warning is variable read length, the expected signature of trimming.
fastqc data/ecoli_R1.trimmed.fastq.gz data/ecoli_R2.trimmed.fastq.gz -o qc/


# ------------------------------------------------------------------------------
# 4. Index the reference genome
# ------------------------------------------------------------------------------
# Build the BWA index so reads can be located quickly along the reference.
bwa index reference/ecoli_ref.fasta


# ------------------------------------------------------------------------------
# 5. Align the trimmed reads to the reference
# ------------------------------------------------------------------------------
# BWA-MEM places each read pair on the reference, producing a SAM file.
bwa mem reference/ecoli_ref.fasta \
        data/ecoli_R1.trimmed.fastq.gz data/ecoli_R2.trimmed.fastq.gz \
        > results/ecoli_aligned.sam


# ------------------------------------------------------------------------------
# 6. Sort and compress the alignment
# ------------------------------------------------------------------------------
# Convert SAM to a coordinate-sorted, compressed BAM (ready for variant calling).
samtools sort results/ecoli_aligned.sam -o results/ecoli_sorted.bam

# Index the sorted BAM for fast positional access.
samtools index results/ecoli_sorted.bam

# Alignment statistics.
samtools flagstat results/ecoli_sorted.bam
# Result: 77.3% of reads mapped. The mapping rate reflects strain divergence
# (O157 sample vs K-12 reference); reads from O157-specific regions have no
# counterpart in the reference and do not map.


# ------------------------------------------------------------------------------
# 7. Call variants
# ------------------------------------------------------------------------------
# mpileup computes the evidence at each position; call decides where variants
# are and writes the compressed VCF.
bcftools mpileup -f reference/ecoli_ref.fasta results/ecoli_sorted.bam \
  | bcftools call -mv -Oz -o results/ecoli_variants.vcf.gz

# Index the VCF.
bcftools index results/ecoli_variants.vcf.gz

# Result: 52,814 raw variants called. This is the unfiltered set; a rigorous
# analysis would filter by quality (QUAL) and depth (DP) to retain the reliable
# subset.
