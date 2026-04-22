#!/usr/bin/env python

# Script to perform indel quantification around gRNA cut site window.

# Load libraries.
import sys
import re
import pandas as pd
import pysam

# Set variables and constants.
cigar_dict = {
        0: 'M',
        1: 'I',
        2: 'D',
        3: 'N',
        4: 'S',
        5: 'H',
        6: 'P',
        7: '=',
        8: 'X',
        9: 'B'
}

# Function to process cigar string.
def process_cigar(cigar, start, stop):
    cigar_full = []
    for i in cigar:
        cigar_full.append(cigar_dict[i[0]] * i[1])
    out = ''.join(cigar_full)
    return out[start:stop]

# Function to locate cut site given reference and protospacer sequences.
def cut_site(reference, protospacer, window_size, cas = 'cas9'):
    if protospacer not in reference:
        raise Exception("Protospacer sequence not found in reference sequence.")
    if cas == 'cas9':
        offset = 3
    pos = re.search(protospacer, reference).span()
    cut_site = pos[1] - offset
    return cut_site - window_size, cut_site + window_size

# Load Amplicon-seq data in bam format.
bam_file = sys.argv[1]
sample_id = sys.argv[2]
reference = sys.argv[3]
protospacer = sys.argv[4]
bam = pysam.AlignmentFile(bam_file, "rb")
mapped_reads = 0
window_size = 8
start, stop = cut_site(reference, protospacer, window_size, 'cas9')
alleles = []

# Loop through reads and extract alleles within quantification window.
for read in bam.fetch():
    read_start = read.reference_start
    read_end = read_start + read.query_length
    cigar_string = read.cigarstring
    cigar_tuple = read.cigartuples
    sequence = read.query_alignment_sequence
    allele = process_cigar(cigar_tuple, start, stop)
    if read.is_mapped and read.mapping_quality >= 10 and read_start <= start and read_end >= stop and 'S' not in allele:
        mapped_reads += 1
        alleles.append(allele)

alleles = pd.Series(alleles)
alleles_freq = alleles.value_counts()
allele_df = alleles_freq.reset_index()
allele_df.columns = ['allele', 'count']
allele_df['percentage'] = (allele_df['count'] / mapped_reads) * 100
allele_df.to_csv(sample_id + '.csv', index = False)
