from functionsmodule import *
from scipy.stats import binom

sample = sys.argv[1]
modus = sys.argv[2]
path = './output/'

#Load PCR error rate:
mismatch_rate = pd.read_csv(path + sample + '/alignmentsummarymetrics.txt', skiprows=4, sep='\t', header=1)['PF_MISMATCH_RATE'][2].astype(float)

#Load annotated variants
annotated = pd.read_excel(path + sample + '/annotated_SNVs.xlsx', sheet_name='Variant', skiprows=1, engine='openpyxl')
annotated['VAF'] = annotated['Samples'].str.split('AF:', 1).str[1].str.split('_').str[0].astype(float)
annotated['DP'] = annotated['Samples'].str.split('DP:', 1).str[1].str.split(';').str[0].astype(int)
annotated['MDP'] = (annotated['DP'] * annotated['VAF']).astype(float)

#CALCULATE ERROR RATE!!!
annotated['ErrorRate'] = (1 - binom.cdf(k=annotated['MDP'], n=annotated['DP'], p=mismatch_rate)).astype(float).round(6)
#!!!!!!!!!!!!!!!!!!!!

#Initiate filter counts df
filtereddf = pd.DataFrame(columns=['filter', 'variants_remaining'])
filtereddf.loc[len(filtereddf)] = ['start', len(annotated)]


###START FILTERING & recording filtered variants
#Synonymous
filteredvariants1 = annotated[annotated['Sequence Ontology'] == 'synonymous_variant'].copy() #Save synonymous variants in other df
filteredvariants1['FilterReason'] = 'Synonymous'
annotated = annotated[annotated['Sequence Ontology'] != 'synonymous_variant'] #Filter synonymous
filtereddf.loc[len(filtereddf)] = ['Synonymous', len(annotated)]

#3' variants
filteredvariants2 = annotated[annotated['Sequence Ontology'] == '3_prime_UTR_variant'].copy()
filteredvariants2['FilterReason'] = "3' UTR"
annotated = annotated[annotated['Sequence Ontology'] != '3_prime_UTR_variant'] #Filter 3' UTR variants
filtereddf.loc[len(filtereddf)] = ["3' UTR", len(annotated)]

#5' variants
filteredvariants3 = annotated[annotated['Sequence Ontology'] == '5_prime_UTR_variant'].copy()
filteredvariants3['FilterReason'] = "5' UTR"
annotated = annotated[annotated['Sequence Ontology'] != '5_prime_UTR_variant'] #Filter 5' UTR variants
filtereddf.loc[len(filtereddf)] = ["5' UTR", len(annotated)]

#Intronic
filteredvariants4 = annotated[annotated['Sequence Ontology'] == 'intron_variant'].copy()
filteredvariants4['FilterReason'] = "Intronic"
annotated = annotated[annotated['Sequence Ontology'] != 'intron_variant'] #Filter Intronic variants
filtereddf.loc[len(filtereddf)] = ['Intronic', len(annotated)]

#2kb upstream
filteredvariants5 = annotated[annotated['Sequence Ontology'] == '2kb_upstream_variant'].copy()
filteredvariants5['FilterReason'] = "2kb upstream"
annotated = annotated[annotated['Sequence Ontology'] != '2kb_upstream_variant'] #Filter 2kb Upstream variants
filtereddf.loc[len(filtereddf)] = ['2kb upstream', len(annotated)]

#2kb downstream
filteredvariants6 = annotated[annotated['Sequence Ontology'] == '2kb_downstream_variant'].copy()
filteredvariants6['FilterReason'] = "2kb downstream"
annotated = annotated[annotated['Sequence Ontology'] != '2kb_downstream_variant'] #Filter 2kb Downstream variants
filtereddf.loc[len(filtereddf)] = ['2kb downstream', len(annotated)]

#PoN
if modus == 'sample':
    pass
elif modus == 'PoN':
    annotated['mut'] = annotated['Chrom'].str[3:] + '_' + annotated['Pos'].astype(str)+'_' + annotated['Reference allele'] +'>' + annotated['Alternate allele']
    dfPON = pd.read_csv('./PoN/BLACKLIST.txt', sep='\t', header=None) #Load blacklist 
    dfPON['mut'] = dfPON[0].astype(str) + '_' + dfPON[1].astype(str)+'_' + dfPON[2] +'>' + dfPON[3]
    
    filteredvariants7 = annotated[annotated.mut.isin(dfPON['mut'])].drop(['mut'], axis=1) .copy() #Save blacklisted variants in other df
    filteredvariants7['FilterReason'] = 'PoN'
    annotated = annotated[~annotated.mut.isin(dfPON['mut'])].drop(['mut'], axis=1) #Filter out blacklisted variants
    
    filtereddf.loc[len(filtereddf)] = ['PoN', len(annotated)]
else:
    print("Wrong input")

#VAF < 0.5%
filteredvariants8 = annotated[annotated['VAF'] < 0.005].copy()
filteredvariants8['FilterReason'] = "VAF < 0.5%"
annotated = annotated[annotated['VAF'] >= 0.005]
filtereddf.loc[len(filtereddf)] = ['VAF < 0.5%', len(annotated)]    
    
#VAF < 1.0%
filteredvariants9 = annotated[annotated['VAF'] < 0.01].copy()
filteredvariants9['FilterReason'] = "VAF < 1%"
annotated = annotated[annotated['VAF'] >= 0.01]
filtereddf.loc[len(filtereddf)] = ['VAF < 1%', len(annotated)]

#in SNP db > 1% AF
filteredvariants10 = annotated.loc[(annotated['AF'].fillna(0) > 0.01) | (annotated['EUR AF'].fillna(0) > 0.01) |
                            (annotated['Global AF'].fillna(0) > 0.01) | (annotated['Non-Fin Eur AF'].fillna(0) > 0.01)].copy()
filteredvariants10['FilterReason'] = "In SNP db > 1%"
annotated = annotated.loc[(annotated['AF'].fillna(0) <= 0.01) | (annotated['EUR AF'].fillna(0) <= 0.01) |
                            (annotated['Global AF'].fillna(0) <= 0.01) | (annotated['Non-Fin Eur AF'].fillna(0) <= 0.01)] #Filter out SNPS with db VAF of > 1%
filtereddf.loc[len(filtereddf)] = ['>1% AF in SNP DB', len(annotated)]   
 
#MNVs:
filteredvariants11 = annotated[annotated['Reference allele'].str.len() | annotated['Alternate allele'].str.len() > 5].copy()
filteredvariants11['FilterReason'] = "MNV (> 5bp)"
annotated = annotated[annotated['Reference allele'].str.len() | annotated['Alternate allele'].str.len() <= 5]
filtereddf.loc[len(filtereddf)] = ["MNV (> 5bp)", len(annotated)]

  
if modus == 'sample':
    filteredvariantslist = [filteredvariants1, filteredvariants2, filteredvariants3, 
                            filteredvariants4, filteredvariants5, filteredvariants6,
                            filteredvariants8, filteredvariants9, filteredvariants10, filteredvariants11]    
elif modus == 'PoN':    
    filteredvariantslist = [filteredvariants1, filteredvariants2, filteredvariants3, 
                            filteredvariants4, filteredvariants5, filteredvariants6, filteredvariants7,
                            filteredvariants8, filteredvariants9, filteredvariants10, filteredvariants11]    

filteredvariants = pd.concat(filteredvariantslist)   




### Saving filtered variants, removed variants and filter counts files
if modus == 'PoN':
    filtereddf.to_csv(path + sample + '/filtering_info_PoN.txt', sep='\t', index = False, header= True)
    annotated.to_excel(path + sample + '/filtered_variants_PoN.xlsx', index=False)
    filteredvariants.to_excel(path + sample + '/removed_variants_PoN.xlsx', index=False)
    
elif modus == 'sample':
    filtereddf.to_csv(path + sample + '/filtering_info.txt', sep='\t', index = False, header= True)
    annotated.to_excel(path + sample + '/filtered_variants.xlsx', index=False)
    filteredvariants.to_excel(path + sample + '/removed_variants.xlsx', index=False)
else: 
    print('wrong input for post filtering...')
    

