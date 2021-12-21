"""
This module automatically converts sorted vardict output to standard .vcf format.

Author: Nick Veltmaat
Date: 19-11-2021
"""

import pandas as pd
import glob
import sys

if len(sys.argv) != 3:
  print("Usage:\t" + sys.argv[0] + "\t<input_sorted_vardict_file_path>\t<output_path_filename>")
  exit(0)

file = glob.glob(str(sys.argv[1]))
data = pd.read_csv(file[0], sep='\t', header=None)

df = data.rename(columns={1: '#CHROM', 3: 'POS', 5:"REF", 6:"ALT"})
df[7] = 'DP=' + df[7].astype(str)
df[8] = 'ADP=' + df[8].astype(str)
df[9] = 'RFW=' + df[9].astype(str)
df[10] = 'RREV=' + df[10].astype(str)
df[11] = 'AFW=' + df[11].astype(str)
df[12] = 'AREV=' + df[12].astype(str)
df[13] = 'GT=' + df[13].astype(str)
df[14] = 'AF=' + df[14].astype(str)
df[15] = 'BIAS=' + df[15].astype(str).str.replace(';','|')
df[16] = 'PMEAN=' + df[16].astype(str)
df[17] = 'PSTD=' + df[17].astype(str)
df[18] = 'QMean=' + df[18].astype(str)
df[19] = 'QStd=' + df[19].astype(str)
df[20] = 'MQ=' + df[20].astype(str)
df[21] = 'SN=' + df[21].astype(str)
df[22] = 'HiAF=' + df[22].astype(str)
df[23] = 'ExAF=' + df[23].astype(str)
df[24] = 'SHIFT3=' + df[24].astype(str)
df[25] = 'MSI=' + df[25].astype(str)
df[26] = 'MSINT=' + df[26].astype(str)
df[27] = 'NM=' + df[27].astype(str)
df[28] = 'HiCnt=' + df[28].astype(str)
df[29] = 'HiCov=' + df[29].astype(str)
df[30] = '5pFS=' + df[30].astype(str)
df[31] = '3pFS=' + df[31].astype(str)
df[32] = 'Seg=' + df[32].astype(str)
df[33] = 'VarType=' + df[33].astype(str)

df = df.drop([0, 2, 4, 34, 35], axis=1)

df['INFO'] = df[df.columns[4:]].apply(    lambda x: ';'.join(x.astype(str)),    axis=1)
df = df.filter(['#CHROM', 'POS', 'REF', 'ALT', 'INFO'])
df['ID'], df['QUAL'], df['FILTER'] = ['.', '.', '.']
df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]

file1 = open(str(sys.argv[2]),"w",newline='')
file1.write('''##fileformat=VCFv4.3\n''')
file1.write('''##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth, total Coverage">\n''')
file1.write('''##INFO=<ID=ADP,Number=1,Type=Integer,Description="No. of reads supporting alternative allele">\n''')
file1.write('''##INFO=<ID=RFW,Number=1,Type=Integer,Description="No. of reads in forward orientation supporting reference allele">\n''')
file1.write('''##INFO=<ID=RREV,Number=1,Type=Integer,Description="No. of reads in reverse orientation supporting reference allele">\n''')
file1.write('''##INFO=<ID=AFW,Number=1,Type=Integer,Description="No. of reads in forward orientation supporting alternative allele">\n''')
file1.write('''##INFO=<ID=AREV,Number=1,Type=Integer,Description="No. of reads in reverse orientation supporting alternative allele">\n''')
file1.write('''##INFO=<ID=GT,Number=1,Type=Flag,Description="Genotype">\n''')
file1.write('''##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency in fraction [0-1]">\n''')
file1.write('''##INFO=<ID=BIAS,Number=1,Type=Flag,Description="Whether there’s strand bias. It consists two numbers. First number is for reference allele. 2nd number is for alternative allele. Numbers are 0-2. 0 means not enough to determine bias. 1 means only one orientation observed. 2 means both orientations observed. 2:1 would indicate strand bias, but not 1:1.">\n''')
file1.write('''##INFO=<ID=PMEAN,Number=1,Type=Float,Description="The mean position in reads of alternative allele. Smaller number would suggest false positives">\n''')
file1.write('''##INFO=<ID=PSTD,Number=1,Type=Integer,Description="Indicate whether the position in reads are the same. 0 means they are the same, 1 mean they are different">\n''')
file1.write('''##INFO=<ID=QMean,Number=1,Type=Float,Description="The mean base quality for alternative allele">\n''')
file1.write('''##INFO=<ID=QStd,Number=1,Type=Integer,Description="Indicate whether the base quality in reads are the same. 0 means they are the same, 1 mean they are different">\n''')
file1.write('''##INFO=<ID=MQ,Number=1,Type=Float,Description="The mean mapping quality for reads supporting alternative allele">\n''')
file1.write('''##INFO=<ID=SN,Number=1,Type=Flag,Description="Signal to noise ratio. The higher the number (>1.5), the more reliable the calls">\n''')
file1.write('''##INFO=<ID=HiAF,Number=1,Type=Float,Description="Allele frequency if only high base quality reads are used">\n''')
file1.write('''##INFO=<ID=ExAF,Number=1,Type=Float,Description="Extra allele frequency recovered from local realignment">\n''')
file1.write('''##INFO=<ID=SHIFT3,Number=1,Type=Integer,Description="No. of bases the Indel can be shifted 3’ with equivalent alignment">\n''')
file1.write('''##INFO=<ID=MSI,Number=1,Type=Integer,Description="Whether there’s microsatellite in sequence context">\n''')
file1.write('''##INFO=<ID=MSINT,Number=1,Type=Integer,Description="Number of bp per unit for MSI. 1 would indicate homopolymer, 2 for di-nucleotide repeat, and so on… ">\n''')
file1.write('''##INFO=<ID=NM,Number=1,Type=Float,Description="Mean number of mismatches in the reads (excluding indels) supporting alternative allele">\n''')
file1.write('''##INFO=<ID=HiCnt,Number=1,Type=Integer,Description="No. of reads with high base quality">\n''')
file1.write('''##INFO=<ID=HiCov,Number=1,Type=Integer,Description="No. of coverage with high base quality">\n''')
file1.write('''##INFO=<ID=5pFS,Number=1,Type=Flag,Description="20bp flanking the variants at 5’">\n''')
file1.write('''##INFO=<ID=3pFS,Number=1,Type=Flag,Description="20bp flanking the variants at 3’">\n''')
file1.write('''##INFO=<ID=Seg,Number=1,Type=Float,Description="The genomic segment variant is called">\n''')
file1.write('''##INFO=<ID=VarType,Number=1,Type=Flag,Description="The type of variants. Values are: SNV, MNV, Insertion, Deletion, and Complex">\n''')

for i in df['#CHROM'].unique():
    file1.write("##contig=<ID="+str(i)+">\n")
# file1.write(df.to_string(index = False, justify='left', colsp))
file1.write(df.to_csv(index = False, sep='\t'))      
file1.close()
