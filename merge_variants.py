from functionsmodule import *

sample = sys.argv[1]
data = sys.argv[2]
path = './output/'
sites = glob.glob(path + sample + '/' + data)

#Load sites.txt
df = pd.read_csv(sites[0], sep='\t', header=None)
df[4] = df[4].astype(str)
df[4] = df[4].apply(lambda x: x.zfill(4))
df['loc'] = df[0].astype(str) + '_' + df[1].astype(str) + '_' + df[2] + '>' + df[3]
df.columns = ["CHROM", "POS", 'REF', 'ALT', "TOOLS",'loc']

#Stripped Mutect2 output:
dfm2 = read_vcf(path + sample +'/0001.vcf')
dfm2 = vcf_info_mutect(dfm2).rename(columns = {'AF':'AFm2', 'RD':'RDm2', 'MDP':'MDPm2'})
dfm2['loc'] = dfm2['CHROM'] + '_' + dfm2['POS'].astype(str) + '_' + dfm2['REF'] + '>' + dfm2['ALT']
dfm2 = dfm2.drop(columns =['CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL',  'INFO', 'FORMAT', 'FILTER']).drop(dfm2.columns[9], axis=1)

#Stripped LoFreq output
dflf = read_vcf(path + sample +'/0002.vcf')
dflf = vcf_info_lofreq(dflf).rename(columns = {'AF':'AFlf', 'RD':'RDlf', 'MDP':'MDPlf'})
dflf['loc'] = dflf['CHROM'] + '_' + dflf['POS'].astype(str) + '_' + dflf['REF'] + '>' + dflf['ALT']
dflf = dflf.drop(columns =['CHROM', 'POS','REF', 'ALT', 'ID', 'FILTER', 'INFO', 'QUAL'])

#Stripped VarDict output
dfvd = read_vcf(path+ sample +'/0003.vcf')
dfvd = vcf_info_vardict(dfvd).rename(columns = {'AF':'AFvd', 'RD':'RDvd', 'MDP':'MDPvd'})
dfvd['loc'] = dfvd['CHROM'] + '_' + dfvd['POS'].astype(str) + '_' + dfvd['REF'] + '>' + dfvd['ALT']
dfvd = dfvd.drop(columns =['CHROM', 'POS', 'REF', 'ALT', 'ID', 'FILTER', 'QUAL', 'INFO'])

#Stripped SiNVICT output
dfsv = read_vcf(path+ sample +'/0000.vcf')
dfsv = vcf_info_sinvict(dfsv).rename(columns = {'AF':'AFsv', 'RD':'RDsv', 'MDP':'MDPsv'})
dfsv['loc'] = dfsv['CHROM'] + '_' + dfsv['POS'].astype(str) + '_' + dfsv['REF'] + '>' + dfsv['ALT']
dfsv = dfsv.drop(columns =['CHROM', 'POS', 'REF', 'ALT', 'ID', 'FILTER', 'QUAL', 'INFO', 'FORMAT', 'calls_level1_sorted.sinvict'])

#Merging dfs per tool
df_merged = pd.merge(pd.merge(pd.merge(pd.merge(dfm2, dflf, on='loc', how='outer'), 
                                       dfsv, on='loc', how='outer'), 
                              dfvd, on='loc', how='outer'), 
                     df, on='loc', how='outer')
df_merged["AF_mean"] = df_merged[["AFm2", "AFlf", "AFsv", "AFvd"]].apply(pd.to_numeric, args=['coerce']).mean(axis=1, skipna = True).round(4)
df_merged["DP_mean"] = df_merged[["RDm2", "RDlf", "RDsv", "RDvd"]].apply(pd.to_numeric, args=['coerce']).mean(axis=1, skipna = True).round(0).astype(int)
df_merged = df_merged[['CHROM', 'POS','REF','ALT','TOOLS','AF_mean','DP_mean']] #drop columns
df_merged['info'] = 'SV-M2-LF-VD:' + df_merged['TOOLS'] + '_AF:' + df_merged['AF_mean'].astype(str) + '_DP:' + df_merged['DP_mean'].astype(str)
df_merged


#Write vcfs_merged.txt
if data == 'sites.txt':
    df_merged.drop(columns =['AF_mean', 'DP_mean', 'TOOLS']).to_csv(path + sample + '/vcfs_merged.txt', sep='\t', index = False, header= False)

#Remove this part: everything will be filtered at the end
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
