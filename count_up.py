#!/usr/bin/env python3

# this is a python script template
# this next line will download the file using curl

gff="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
fasta="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

import os, gzip, itertools, csv, re

# this is code which will parse FASTA files
# define what a header looks like in FASTA format

def isheader(line):
    return line[0] == '>'


def aspairs(f):
    seq_id = ''
    sequence = ''
    for header, group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence


# download the file

if not os.path.exists(gff):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")
    
with gzip.open(gff,"rt") as fh:
    # now add code to process this
    gff = csv.reader(fh,delimiter="\t")
    for row in gff:
        if row[0].startswith("#"):
            continue
        print(row[3],row[6])



# count up and print the number of genes in this file
# compute the total length of thegenes


n = 0 
Total_length = 0 
gff2 = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
with gzip.open(gff2, "rt") as fh:
  gff = csv.reader(fh, delimiter="\t")
  for row in gff:
    if row[0].startswith("Chromosome") and row[2]=="gene":
      Total_length += int(row[4])-int(row[3])

if row[0].startswith("###"):
  n += 1

  print("The total number of genes is: {}".format(n))
  print("The total length of the genes is: {}".format(Total_length))


import itertools
import gzip
import sys 
import re


# define what a header looks like in FASTA format

def isheader(line):
    return line[0] == '>'


def aspairs(f):
    seq_id = ''
    sequence = ''
    for header, group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence


# to compute the total length of the genome 

with gzip.open(fasta, "rt") as fh: 
# with open(filename, "r") as f:


  pairs = aspairs(fh)
  seqs = dict(pairs)


  # to iterate through the sequences

  n = 0
  print("The total length of the genome is: {}".format(len(seqs["Chromosome"])))

  # print out the % of genome which is coding
print("The percentage of the E.coli K12 genome which is coding: {:.2f}%".format(Total_length/len(seqs["Chromosome"])*100))

  
