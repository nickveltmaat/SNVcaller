"""
This module loads 'sites.txt' from a specific sample in the output folder,
Generates a Venn diagram of amount of mutations called per tool,
Generates Histogram of amount of mutations with a certain Read Depth
Generates Histogram of amount of mutations with a certain VAF (Variant Allele Frequency)

Author: Nick Veltmaat
Date: 17-11-2021
"""

import io
import os
import re
import sys
import math
import glob
import pandas_bokeh
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from functools import reduce
from venn import venn

## Load data (sites.txt)
sample = sys.argv[1]
data = sys.argv[2]
path = './output/'

if data == 'sites.txt':
    pontype = ""
elif data == 'sites_PoN.txt':
    pontype = "PoN"
else:
  print("Wrong input")

sites = glob.glob(path + sample +'/' + data)
df = pd.read_csv(sites[0], sep='\t', header=None)
df[4] = df[4].astype(str)
df[4] = df[4].apply(lambda x: x.zfill(4))
df['loc'] = df[0] + '_' + df[1].astype(str) + '_' + df[2] + '>' + df[3]
df.columns = ["CHROM", "POS", 'REF', 'ALT', "TOOLS",'loc']


## Extracting data from vcfs
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

#Stripped Mutect2 output:
dfm2 = read_vcf(path + sample +'/0001.vcf')
dfm2 = dfm2[dfm2['FILTER']=='PASS']
dfm2['DPm2'] = dfm2[dfm2.columns[9]].str.split(':', 1).str[1].str.strip().str.split(':', 1).str[1].str.strip().str.split(':', 1).str[1].str.split(':', 0).str[0].astype(int)
dfm2['AFm2'] = dfm2[dfm2.columns[9]].str.split(':', 1).str[1].str.strip().str.split(':', 1).str[1].str.strip().str.split(':', 1).str[0].str.strip().astype(float)
dfm2['loc'] = dfm2['CHROM'] + '_' + dfm2['POS'].astype(str) + '_' + dfm2['REF'] + '>' + dfm2['ALT']
dfm2 = dfm2.drop(columns =['CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
dfm2 = dfm2.drop(dfm2.columns[0], axis=1)

#Stripped LoFreq output
dflf = read_vcf(path + sample +'/0002.vcf')
dflf['DPlf'] = dflf['INFO'].str.split('DP=', 1).str[1].str.strip().str.split(';', 1).str[0].str.strip().astype(int)
dflf['AFlf'] = dflf['INFO'].str.split('AF=', 1).str[1].str.strip().str.split(';', 1).str[0].str.strip().astype(float)
dflf['loc'] = dflf['CHROM'] + '_' + dflf['POS'].astype(str) + '_' + dflf['REF'] + '>' + dflf['ALT']
dflf = dflf.drop(columns =['CHROM', 'POS','REF', 'ALT', 'ID', 'FILTER', 'INFO', 'QUAL'])

#Stripped VarDict output
dfvd = read_vcf(path+ sample +'/0003.vcf')
dfvd['DPvd'] = dfvd['INFO'].str.split('DP=', 1).str[1].str.strip().str.split(';', 1).str[0].str.strip().astype(int)
dfvd = dfvd.loc[dfvd['DPvd'] > int(sys.argv[3]) ] # Filter on RD (since VarDict doesn't have a parameter for this)
dfvd['AFvd'] = dfvd['INFO'].str.split('AF=', 1).str[1].str.strip().str.split(';', 1).str[0].str.strip().astype(float)
dfvd['loc'] = dfvd['CHROM'] + '_' + dfvd['POS'].astype(str) + '_' + dfvd['REF'] + '>' + dfvd['ALT']
dfvd = dfvd.drop(columns =['CHROM', 'POS', 'REF', 'ALT', 'ID', 'FILTER', 'QUAL', 'INFO'])

#Stripped SiNVICT output
dfsv = read_vcf(path+ sample +'/0000.vcf')
dfsv['DPsv'] = dfsv['calls_level1_sorted.sinvict'].str.split(':', 1).str[0]
dfsv['AFsv'] = dfsv['calls_level1_sorted.sinvict'].str.split(':', 1).str[1]
dfsv['loc'] = dfsv['CHROM'] + '_' + dfsv['POS'].astype(str) + '_' + dfsv['REF'] + '>' + dfsv['ALT']
dfsv = dfsv.drop(columns =['CHROM', 'POS', 'REF', 'ALT', 'ID', 'FILTER', 'QUAL', 'INFO', 'FORMAT', 'calls_level1_sorted.sinvict'])

dataframes = [df, dfm2, dflf, dfvd, dfsv]

df_merged = pd.merge(pd.merge(pd.merge(pd.merge(dfm2, dflf, on='loc', how='outer'), 
                                       dfsv, on='loc', how='outer'), 
                              dfvd, on='loc', how='outer'), 
                     df, on='loc', how='outer')
df_merged["AF_mean"] = df_merged[["AFm2", "AFlf", "AFsv", "AFvd"]].apply(pd.to_numeric, args=['coerce']).mean(axis=1, skipna = True).round(4)
df_merged["DP_mean"] = df_merged[["DPm2", "DPlf", "DPsv", "DPvd"]].apply(pd.to_numeric, args=['coerce']).mean(axis=1, skipna = True).round()
df_merged = df_merged.drop(columns =['loc', 'DPm2', 'AFm2', 'DPlf', 'AFlf', 'DPsv', 'AFsv', 'DPvd', 'AFvd']).dropna()
df_merged['POS'] = df_merged['POS'].astype(int)
df_merged['DP_mean'] = df_merged['DP_mean'].astype(int)
df_merged['info'] = 'SV-M2-LF-VD:' + df_merged['TOOLS'] + '_AF:' + df_merged['AF_mean'].astype(str) + '_DP:' + df_merged['DP_mean'].astype(str) 

#Export "vcfs_merged.txt"
if data == 'sites.txt':
  df_merged.drop(columns =['AF_mean', 'DP_mean', 'TOOLS']).to_csv(path + sample + '/vcfs_merged.txt', sep='\t', index = False, header= False)
elif data == 'sites_PoN.txt':
  df_merged = df_merged[df_merged['REF'].notna()]
  df_merged['mut'] = df_merged['CHROM'] + '_' + df_merged['POS'].astype(str)+'_' + df_merged['REF'] +'>' + df_merged['ALT']
  dfPON = pd.read_csv('./PoN/BLACKLIST.txt', sep='\t', header=None)
  dfPON['mut'] = dfPON[0].astype(str) + '_' + dfPON[1].astype(str)+'_' + dfPON[2] +'>' + dfPON[3]
  df_merged = df_merged[~df_merged.mut.isin(dfPON['mut'])].drop(['mut'], axis=1)
  df_merged.drop(columns =['AF_mean', 'DP_mean', 'TOOLS']).to_csv(path + sample + '/vcfs_merged_PoN.txt', sep='\t', index = False, header= True)
else:
  print("Wrong input")
  
#Venn diagram
dfvenn = df_merged.copy()

dfvenn['SV'] = dfvenn['TOOLS'].str[0:1]
dfvenn['M2'] = dfvenn['TOOLS'].str[1:2]
dfvenn['LF'] = dfvenn['TOOLS'].str[2:3]
dfvenn['VD'] = dfvenn['TOOLS'].str[3:4]
dfvenn['loc'] = dfvenn['CHROM'] + '_' + dfvenn['POS'].astype(str) + '_' + dfvenn['REF'] + '>' + dfvenn['ALT']

dfSV = dfvenn.loc[dfvenn['SV'] == '1']
SVlist = list(dfSV['loc'])
dfM2 = dfvenn.loc[dfvenn['M2'] == '1']
M2list = list(dfM2['loc'])
dfLF = dfvenn.loc[dfvenn['LF'] == '1']
LFlist = list(dfLF['loc'])
dfVD = dfvenn.loc[dfvenn['VD'] == '1']
VDlist = list(dfVD['loc'])

SNVcalls = {"SiNVICT": {i for i in SVlist},
            "Mutect2": {i for i in M2list},
            "LoFreq": {i for i in LFlist},
            "VarDict": {i for i in VDlist}}

fig = venn(SNVcalls
     #, fmt="{percentage:.1f}%"
    ).figure
fig.suptitle(sys.argv[1] + "  " + pontype, fontsize=15)
plt.xlabel('Total # of mutations: '+ str(len(dfvenn)), fontsize=12)
fig.savefig(path + sample + '/venn'+ pontype+ '.png')   # save the figure to file
plt.close(fig) 


#Histograms
def roundup(x):
  """
  This function rounds up a number to the first higher 200-fold
  
  Argument: Float or Int
  
  Returns 200-fold int
  """
  return int(math.ceil(x / 200.0)) * 200
    
  
hist_AF = df_merged['AF_mean'].plot_bokeh(
    kind="hist",
    bins=np.linspace(0, 1, 101),
    histogram_type="sidebyside",
    vertical_xlabel=True,
    hovertool=True,
    title="Allele Frequency " + sample,
    ylabel = '# Of mutations',
    xlabel = 'Mean Allele Frequency per mutation, calculated by Mutect2, LoFreq & VarDict',
    legend = None,
    line_color="black",
    show_figure = False)

hist_AF2 = df_merged['AF_mean'].plot_bokeh(
    kind="hist",
    color='darkblue',
    bins=np.linspace(0, 0.01, 21),
    histogram_type="sidebyside",
    vertical_xlabel=True,
    hovertool=True,
    title="Allele Frequency " + sample,
    ylabel = '# Of mutations',
    xlabel = 'Mean Allele Frequency per mutation, calculated by Mutect2, LoFreq & VarDict, Zoomed in on 0 - 0.001',
    legend = None,
    line_color="black",
    show_figure = False)

hist_DP = df_merged['DP_mean'].plot_bokeh(
    kind="hist",
    color='green',
    bins=np.linspace(0, roundup(df_merged['DP_mean'].max()), int(roundup(df_merged['DP_mean'].max())/200+1)),
    histogram_type="sidebyside",
    vertical_xlabel=True,
    hovertool=True,
    title='Read Depth '+ sample,
    ylabel = '# Of mutations',
    xlabel = 'Mean Read Depth per mutation, calculated by Mutect2, LoFreq & VarDict',
    legend = None,
    line_color="black",
    show_figure = False)

pandas_bokeh.output_file(path + sample +"/Mean_Depth_&_AF_per_Mutation_"+ pontype + ".html")  
pandas_bokeh.plot_grid([[hist_AF], [hist_AF2], [hist_DP]], plot_width=1500, plot_height=440)
