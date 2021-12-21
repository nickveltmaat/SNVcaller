"""
This module converts sorted sinvict output to standard .vcf format.
Two arguments required to be given when calling this module:
  #1: path to sorted bcftools isec or sinvict file as input
  #2: Path to location and name for newly generated .vcf file

Author: Nick Veltmaat
Date: 20-12-2021
"""

import sys
import subprocess
if len(sys.argv) != 5:
  print("Usage:\t" + sys.argv[0] + "\t<input_sorted_bcftools_output_file>\t<output_filename>\t<reference_fasta>\t<input_type: 'sinvict' or 'bcftools_isec'>")
  exit(0)
  
outfile = open(sys.argv[2], "w")
outfile.write("##fileformat=VCFv4.3\n")

chroms = []
with open(sys.argv[1]) as infile:
  for line in infile:
    line = line.rstrip()
    tokens = line.split()
    chromosome = tokens[0]
    chroms.append(chromosome)

chroms2 = set(chroms)    
for i in chroms2:
    outfile.write("##contig=<ID="+str(i)+">\n")

outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

with open(sys.argv[1]) as infile:
  for line in infile:
    line = line.rstrip()
    tokens = line.split()
    chromosome = tokens[0]
    position = tokens[1]
    if sys.argv[4] == 'bcftools_isec':
      ref = tokens[2]
      alt = tokens[3]
    elif sys.argv[4] == 'sinvict':
      ref = tokens[3]
      alt = tokens[5]    
    else: 
      print("ERROR in providing correct input type. Choose between 'sinvict' or 'bcftools_isec'")
    if alt[0] == "+" :
      #insertion
      alt = ref + alt[1:]
      outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.upper() + "\t" + alt.upper() + "\t" + "." + "\t" + "PASS" + "\t" + "." + "\n")
    elif alt[0] == "-" :
      #deletion
      position = int(position)
      position -= 1
      fa = subprocess.check_output(["samtools", "faidx", sys.argv[3], chromosome+":"+str(position)+"-"+str(position)])
      base = fa.split("\n".encode())[1]
      base.rstrip()
      ref = base + str(alt[1:]).encode()
      alt = base
      outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.decode("utf-8").upper() + "\t" + alt.decode("utf-8").upper() + "\t" + "." + "\t" + "PASS" + "\t" + "." + "\n")
    else: 
      outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.upper() + "\t" + alt.upper() + "\t" + "." + "\t" + "PASS" + "\t" + "." + "\n")
      
# original copied and modified at 19okt2021 from: https://raw.githubusercontent.com/sfu-compbio/sinvict/master/sinvict_to_vcf.py
