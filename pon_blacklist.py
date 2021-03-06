"""
This module loads 'sites.txt' from a specific sample in the output folder, and loads the previously generated merged_PoN.vcf
merged_PoN will act as a blacklist for sites.txt, essentially filtering out mutations from sites.txt that are present in merged_PoN.vcf

Author: Nick Veltmaat
Date: 1-12-2021
"""
import pandas as pd
import numpy as np
import io
import sys

sample = str(sys.argv[1])
print('sample = '+sample)

def read_vcf(path):
  """
  This function reads a .vcf file and stores it as a df.
  
  Required Argument: path to .vcf file
  
  Returns: Dataframe
  """
  with open(path, 'r') as f:
      lines = [l for l in f if not l.startswith('##')]
  return pd.read_csv(
      io.StringIO(''.join(lines)),
      dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
             'QUAL': str, 'FILTER': str, 'INFO': str},
      sep='\t'
  ).rename(columns={'#CHROM': 'CHROM'})
  
dfPON = pd.read_csv('./PoN/BLACKLIST.txt', sep='\t', header=None)
dfPON['mut'] = dfPON[0].astype(str) + '_' + dfPON[1].astype(str)+'_' + dfPON[2] +'>' + dfPON[3]

df = pd.read_csv('./output/'+ sample +'/sites.txt', sep='\t', header=None)
df[4] = df[4].astype(str)
df[4] = df[4].apply(lambda x: x.zfill(4))
df['mut'] = df[0] + '_' + df[1].astype(str)+'_' + df[2] +'>' + df[3]

df[~df.mut.isin(dfPON['mut'])].drop(['mut'], axis=1).to_csv('./output/'+ sample + '/sites_PoN.txt', sep='\t', index = False, header= False)
