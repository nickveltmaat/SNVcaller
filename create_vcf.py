"""
This module converts sorted sinvict output to standard .vcf format.
Two arguments required to be given when calling this module:
  #1: path to sorted bcftools isec or sinvict file as input
  #2: Path to location and name for newly generated .vcf file

Author: Nick Veltmaat
Date: 20-12-2021
"""
###############
# Include RD & VAF tokens for bcftools isec input

import sys
import subprocess
import os
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

if sys.argv[4] == 'bcftools_isec':
  pass
elif sys.argv[4] == 'sinvict':
  outfile.write('''##FORMAT=<ID=DP,Number=1,Type=String,Description="Read Depth">\n''')
  outfile.write('''##FORMAT=<ID=VAF,Number=1,Type=String,Description="Variant Allele Frequency">\n''')
  outfile.write('''##FORMAT=<ID=MDP,Number=1,Type=String,Description="Number of reads supporting the mutation">\n''')

chroms2 = set(chroms)    
for i in chroms2:
    outfile.write("##contig=<ID="+str(i)+">\n")

outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+os.path.basename(sys.argv[1])+"\n")

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
      vaf = round(float(tokens[7])/100, 3)
      depth = int(tokens[4])
      mutdp = int(tokens[6])
    else: 
      print("ERROR in providing correct input type. Choose between 'sinvict' or 'bcftools_isec'")
    if alt[0] == "+" :
      #insertion
      alt = ref + alt[1:]
      if sys.argv[4] == 'bcftools_isec':
        outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.upper() + "\t" + alt.upper() 
                      + "\t" + "." + "\t" + "PASS" + "\t" + "." + "\t" + "." + "\t" +'.' + "\n")
      elif sys.argv[4] == 'sinvict':
        outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.upper() + "\t" + alt.upper() 
                      + "\t" + "." + "\t" + "PASS" + "\t" + "." + "\t" + "DP:VAF:MDP" + "\t" + str(depth)+':'+str(vaf)+':'+str(mutdp)  + "\n")
    elif alt[0] == "-" :
      #deletion
      position = int(position)
      position -= 1
      fa = subprocess.check_output(["samtools", "faidx", sys.argv[3], chromosome+":"+str(position)+"-"+str(position)])
      base = fa.split("\n".encode())[1]
      base.rstrip()
      ref = base + str(alt[1:]).encode()
      alt = base
      if sys.argv[4] == 'bcftools_isec':
        outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.decode("utf-8").upper() + "\t" + alt.decode("utf-8").upper() 
                      + "\t" + "." + "\t" + "PASS"  + "\t" + "." + "\t" + "." + "\t" + '.' + "\n")
      elif sys.argv[4] == 'sinvict':
        outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.decode("utf-8").upper() + "\t" + alt.decode("utf-8").upper() 
                      + "\t" + "." + "\t" + "PASS"  + "\t" + "." + "\t" + "DP:VAF:MDP" + "\t" + str(depth)+':'+str(vaf)+':'+str(mutdp) + "\n")
    else:
      if sys.argv[4] == 'bcftools_isec':
        outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.upper() + "\t" + alt.upper() + 
        "\t" + "." + "\t" + "PASS"  + "\t" + "." + "\t" + "." + "\t" + '.' + "\n")
      elif sys.argv[4] == 'sinvict':
        outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.upper() + "\t" + alt.upper() + 
        "\t" + "." + "\t" + "PASS"  + "\t" + "." + "\t" + "DP:VAF:MDP" + "\t" + str(depth)+':'+str(vaf)+':'+str(mutdp) + "\n")
      
# original copied and modified at 19okt2021 from: https://raw.githubusercontent.com/sfu-compbio/sinvict/master/sinvict_to_vcf.py
