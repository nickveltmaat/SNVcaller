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

## Venn Diagram
sample = sys.argv[1]
print('sample = '+sample)
sites = glob.glob('./output/'+ sample +'/sites.txt')

df = pd.read_csv(sites[0], sep='\t', header=None)
df[4] = df[4].astype(str)
df[4] = df[4].apply(lambda x: x.zfill(4))

df['SV'] = df[4].str[0:1]
df['M2'] = df[4].str[1:2]
df['LF'] = df[4].str[2:3]
df['VD'] = df[4].str[3:4]
df['loc'] = df[0] + '_' + df[1].astype(str)

dfSV = df.loc[df['SV'] == '1']
SVlist = list(dfSV['loc'])
dfM2 = df.loc[df['M2'] == '1']
M2list = list(dfM2['loc'])
dfLF = df.loc[df['LF'] == '1']
LFlist = list(dfLF['loc'])
dfVD = df.loc[df['VD'] == '1']
VDlist = list(dfVD['loc'])

df.columns = ["CHROM", "POS", 'REF', 'ALT', "TOOLS", 'SV','M2','LF','VD','loc']

SNVcalls = {
    "SiNVICT": {i for i in SVlist},
    "Mutect2": {i for i in M2list},
    "LoFreq": {i for i in LFlist},
    "VarDict": {i for i in VDlist}
}

from venn import venn

fig = venn(SNVcalls
     #, fmt="{percentage:.1f}%"
    ).figure
fig.suptitle(sys.argv[1], fontsize=15)
plt.xlabel('Total # of mutations: '+ str(len(df)), fontsize=12)
fig.savefig('./output/'+ sample +'/venn.png')   # save the figure to file
plt.close(fig) 
len(df)


## Histograms
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
dfm2 = read_vcf('./output/'+ sample +'/0001.vcf')
dfm2 = dfm2[dfm2['FILTER']=='PASS']
dfm2['DP'] = dfm2['INFO'].str.split('DP=', 1).str[1].str.strip().str.split(';', 1).str[0].str.strip().astype(int)
dfm2['AF'] = dfm2[dfm2.columns[9]].str.split(':', 1).str[1].str.strip().str.split(':', 1).str[1].str.strip().str.split(':', 1).str[0].str.strip().astype(float)
dfm2 = dfm2.drop(columns =['REF', 'ALT', 'ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
dfm2 = dfm2.drop(dfm2.columns[2], axis=1)

#Stripped LoFreq output
dflf = read_vcf('./output/'+ sample +'/0002.vcf')
dflf['DP'] = dflf['INFO'].str.split('DP=', 1).str[1].str.strip().str.split(';', 1).str[0].str.strip().astype(int)
dflf['AF'] = dflf['INFO'].str.split('AF=', 1).str[1].str.strip().str.split(';', 1).str[0].str.strip().astype(float)
dflf = dflf.drop(columns =['REF', 'ALT', 'ID', 'FILTER', 'INFO'])

#Stripped VarDict output
dfvd = read_vcf('./output/'+ sample +'/0003.vcf')
dfvd['DP'] = dfvd['INFO'].str.split('DP=', 1).str[1].str.strip().str.split(';', 1).str[0].str.strip().astype(int)
dfvd['AF'] = dfvd['INFO'].str.split('AF=', 1).str[1].str.strip().str.split(';', 1).str[0].str.strip().astype(float)
dfvd = dfvd.drop(columns =['REF', 'ALT', 'ID', 'FILTER', 'QUAL', 'INFO'])

dataframes = [df, dfm2, dflf, dfvd]

df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['CHROM', 'POS'], how='outer'), dataframes)
df_merged["AF_mean"] = df_merged[["AF", "AF_x", "AF_y"]].mean(axis=1)
df_merged["DP_mean"] = df_merged[["DP", "DP_x", "DP_y"]].mean(axis=1).round()
df_merged = df_merged.drop(columns =['TOOLS', 'DP_x', 'AF_x', 'DP_y', 'AF_y', 'DP', 'AF'])

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


hist_DP = df_merged['DP_mean'].plot_bokeh(
    kind="hist",
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

pandas_bokeh.output_file('./output/'+ sample +"/Mean_Depth_&_AF_per_Mutation.html")  
pandas_bokeh.plot_grid([[hist_AF], [hist_DP]], plot_width=1550, plot_height=440)
