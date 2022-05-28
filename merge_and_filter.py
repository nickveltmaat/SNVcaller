import io
import os
import sys
import glob
import numpy as np
import pandas as pd

sample = sys.argv[1]
data = sys.argv[2]
path = './output/'
sites = glob.glob(path + sample + '/' + data)

filtereddf = pd.DataFrame(columns=['filter', 'variants_remaining'])

df = pd.read_csv(sites[0], sep='\t', header=None)
df[4] = df[4].astype(str)
df[4] = df[4].apply(lambda x: x.zfill(4))
df['loc'] = df[0] + '_' + df[1].astype(str) + '_' + df[2] + '>' + df[3]
df.columns = ["CHROM", "POS", 'REF', 'ALT', "TOOLS",'loc']

if data == 'sites.txt':
    pontype = ""
    filtereddf.loc[len(filtereddf)] = ['start', len(df)]
elif data == 'sites_PoN.txt':
    pontype = "PoN"
    nopondf = pd.read_csv(glob.glob(path + sample +'/' + 'sites.txt')[0], sep='\t', header=None)
    filtereddf.loc[len(filtereddf)] = ['start', len(nopondf)]
    filtereddf.loc[len(filtereddf)] = ['PoN', len(df)]
else:
  print("Wrong input")

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
# dfm2 = dfm2[dfm2['FILTER']=='PASS']
dfm2['DPm2'] = dfm2[dfm2.columns[9]].str.split(':', 1).str[1].str.strip().str.split(':', 1).str[1].str.strip().str.split(':', 1).str[1].str.split(':', 0).str[0].astype(int)
dfm2['AFm2'] = dfm2[dfm2.columns[9]].str.split(':', 1).str[1].str.strip().str.split(':', 1).str[1].str.strip().str.split(':', 1).str[0].str.strip().astype(float)
dfm2['loc'] = dfm2['CHROM'] + '_' + dfm2['POS'].astype(str) + '_' + dfm2['REF'] + '>' + dfm2['ALT']
dfm2 = dfm2.drop(columns =['CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL',  'INFO', 'FORMAT'])
dfm2 = dfm2.drop(dfm2.columns[1], axis=1)

#Stripped LoFreq output
dflf = read_vcf(path + sample +'/0002.vcf')
dflf['DPlf'] = dflf['INFO'].str.split('DP=', 1).str[1].str.strip().str.split(';', 1).str[0].str.strip().astype(int)
dflf['AFlf'] = dflf['INFO'].str.split('AF=', 1).str[1].str.strip().str.split(';', 1).str[0].str.strip().astype(float)
dflf['loc'] = dflf['CHROM'] + '_' + dflf['POS'].astype(str) + '_' + dflf['REF'] + '>' + dflf['ALT']
dflf = dflf.drop(columns =['CHROM', 'POS','REF', 'ALT', 'ID', 'FILTER', 'INFO', 'QUAL'])

#Stripped VarDict output
dfvd = read_vcf(path+ sample +'/0003.vcf')
dfvd['DPvd'] = dfvd['INFO'].str.split('DP=', 1).str[1].str.strip().str.split(';', 1).str[0].str.strip().astype(int)
# dfvd = dfvd.loc[dfvd['DPvd'] > int(100) ] # Filter on RD (since VarDict doesn't have a parameter for this)
dfvd['AFvd'] = dfvd['INFO'].str.split('AF=', 1).str[1].str.strip().str.split(';', 1).str[0].str.strip().astype(float)
dfvd['loc'] = dfvd['CHROM'] + '_' + dfvd['POS'].astype(str) + '_' + dfvd['REF'] + '>' + dfvd['ALT']
dfvd = dfvd.drop(columns =['CHROM', 'POS', 'REF', 'ALT', 'ID', 'FILTER', 'QUAL', 'INFO'])

#Stripped SiNVICT output
dfsv = read_vcf(path+ sample +'/0000.vcf')
dfsv['DPsv'] = dfsv['calls_level1_sorted.sinvict'].str.split(':', 1).str[0]
dfsv['AFsv'] = dfsv['calls_level1_sorted.sinvict'].str.split(':', 1).str[1]
dfsv['loc'] = dfsv['CHROM'] + '_' + dfsv['POS'].astype(str) + '_' + dfsv['REF'] + '>' + dfsv['ALT']
dfsv = dfsv.drop(columns =['CHROM', 'POS', 'REF', 'ALT', 'ID', 'FILTER', 'QUAL', 'INFO', 'FORMAT', 'calls_level1_sorted.sinvict'])

#Merging dfs per tool
df_merged = pd.merge(pd.merge(pd.merge(pd.merge(dfm2, dflf, on='loc', how='outer'), 
                                       dfsv, on='loc', how='outer'), 
                              dfvd, on='loc', how='outer'), 
                     df, on='loc', how='outer')
df_merged["AF_mean"] = df_merged[["AFm2", "AFlf", "AFsv", "AFvd"]].apply(pd.to_numeric, args=['coerce']).mean(axis=1, skipna = True).round(4)
df_merged["DP_mean"] = df_merged[["DPm2", "DPlf", "DPsv", "DPvd"]].apply(pd.to_numeric, args=['coerce']).mean(axis=1, skipna = True).round()
df_merged = df_merged.drop(columns =['loc', 'DPm2', 'AFm2', 'DPlf', 'AFlf', 'DPsv', 'AFsv', 'AFvd'])#.dropna()
# df_merged['POS'] = df_merged['POS'].astype(int)
df_merged['DP_mean'] = df_merged['DP_mean'].astype(int)
df_merged['info'] = 'SV-M2-LF-VD:' + df_merged['TOOLS'] + '_AF:' + df_merged['AF_mean'].astype(str) + '_DP:' + df_merged['DP_mean'].astype(str) 
df_merged['FILTER'] = df_merged['FILTER'].fillna('PASS') #Replace NaN with PASS


def filtering(df_merged, filtereddf):
    df_merged = df_merged.loc[(df_merged['FILTER'] == 'PASS')] #Mutect2 PASS filter
    filtereddf.loc[len(filtereddf)] = ['M2 PASS', len(df_merged)]
    df_merged = df_merged[df_merged['DPvd'] > 100] #VarDict RD Filter
    filtereddf.loc[len(filtereddf)] = ['VD RD', len(df_merged)]
    df_merged = df_merged.drop(columns =['FILTER', 'DPvd'])
    return df_merged, filtereddf

if data == 'sites.txt':
    df_merged, filtereddf = filtering(df_merged, filtereddf)

    df_merged.drop(columns =['AF_mean', 'DP_mean', 'TOOLS']).to_csv(path + sample + '/vcfs_merged.txt', sep='\t', index = False, header= False)
    filtereddf.to_csv(path + sample + '/filters.txt', sep='\t', index = False, header= True)

elif data == 'sites_PoN.txt':
    df_merged = df_merged[df_merged['REF'].notna()]    
    df_merged['mut'] = df_merged['CHROM'] + '_' + df_merged['POS'].astype(str)+'_' + df_merged['REF'] +'>' + df_merged['ALT']
    dfPON = pd.read_csv('./PoN/BLACKLIST.txt', sep='\t', header=None)
    dfPON['mut'] = dfPON[0].astype(str) + '_' + dfPON[1].astype(str)+'_' + dfPON[2] +'>' + dfPON[3]
    df_merged = df_merged[~df_merged.mut.isin(dfPON['mut'])].drop(['mut'], axis=1)

    df_merged, filtereddf = filtering(df_merged, filtereddf)
    df_merged['POS'] = df_merged['POS'].astype(int)######
    df_merged.drop(columns =['AF_mean', 'DP_mean', 'TOOLS']).to_csv(path + sample + '/vcfs_merged_PoN.txt', sep='\t', index = False, header= False)
    filtereddf.to_csv(path + sample + '/filters_PoN.txt', sep='\t', index = False, header= True)
else:
    print("Wrong input")
