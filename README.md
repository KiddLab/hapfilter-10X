# hapfilter-10X

This repository contains scripts for filtering candidate mosaic SNVs using 10X data.

## Required Components
Python version v2.x
Python modules: pysam, gzip

## Inputs
A VCF-like file listing candidate SNVs (only first 7 columns are used)
A 10X BAM file (from long ranger) from matched tissue

## Outputs
A crude filter assessment as to whether the 10X data indicates that a site is unlikely to
be real

## Description
The phase information from 10X segregates reads into three classes: haplotype 1, haplotype 2, and
unassigned.  The underlying idea is that, if found at all given depth, real mosaic SNVs should only
occur one  haplotype.   Also, if a mosaic allele makes up a large fraction (>90%) of the reads
at a hapltoype, it is very likely to be a heterozygote. Validations as part of the Brain 
Somatic Mosaicism Network suggest that sites that fail these filters are unlikely to be real.

Given typical read depths, many candidates will have no information in a typical 10X BAM file.


## Usage

'''
python2 hapfilter-10X/annotate-vcf-with-10x.py \
--vcf a-vcf-like-file.vcf \
--bam /path/to/10x/longranger/bam/tissue.bam \
--out a-vcf-like-file.10X-table.txt
'''
If input vcf is gzipped use options --gzip.  To only consider sites that are PASS in the VCF file use
options --pass.  This can greatly speed up execution time.  For large files, consider running in batches.

The output is a text file with columns reporting allele count information for each haplotype.



