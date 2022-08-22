from functionsmodule import *

sample = sys.argv[1]
modus = sys.argv[2]
path = './output/'

#Load annotated variants
annotated = pd.read_excel(path + sample + '/annotated_SNVs.xlsx', sheet_name='Variant', skiprows=1, engine='openpyxl')

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
    filteredvariantslist = [filteredvariants1, filteredvariants2, filteredvariants3, 
                            filteredvariants4, filteredvariants5, filteredvariants6]
elif modus == 'PoN':
    annotated['mut'] = annotated['Chrom'].str[3:] + '_' + annotated['Pos'].astype(str)+'_' + annotated['Reference allele'] +'>' + annotated['Alternate allele']
    dfPON = pd.read_csv('./PoN/BLACKLIST.txt', sep='\t', header=None) #Load blacklist 
    dfPON['mut'] = dfPON[0].astype(str) + '_' + dfPON[1].astype(str)+'_' + dfPON[2] +'>' + dfPON[3]
    
    filteredvariants7 = annotated[annotated.mut.isin(dfPON['mut'])].drop(['mut'], axis=1) .copy() #Save PoN variants in other df
    filteredvariants7['FilterReason'] = 'PoN'
    annotated = annotated[~annotated.mut.isin(dfPON['mut'])].drop(['mut'], axis=1) #Filter out blacklisted mutations
    
    filtereddf.loc[len(filtereddf)] = ['PoN', len(annotated)]
    filteredvariantslist = [filteredvariants1, filteredvariants2, filteredvariants3, 
                            filteredvariants4, filteredvariants5, filteredvariants6, filteredvariants7]
else:
    print("Wrong input")

filteredvariants = pd.concat(filteredvariantslist)


##Counting Vardict unique calls (not filtered out ... yet)
#Vardict unique
#vd = annotated.copy()
#vd['tools'] = vd['Samples'].str.split('VD:', 1).str[1].str.split('_').str[0]
#vd[vd['tools'] != '0001']
#filtereddf.loc[len(filtereddf)] = ['Vardict unique', len(vd[vd['tools'] != '0001'])]

##Counting Variants lower than 0.5%  &  1%
tempdf = annotated.copy()
tempdf['VAF'] = tempdf['Samples'].str.split('AF:', 1).str[1].str.split('_').str[0].astype(float)
tempdf = tempdf[tempdf['VAF'] >= 0.005]
filtereddf.loc[len(filtereddf)] = ['VAF < 0.5%', len(tempdf)]
tempdf = tempdf[tempdf['VAF'] >= 0.01]
filtereddf.loc[len(filtereddf)] = ['VAF < 1%', len(tempdf)]

#Counting Variants in SNP databases with AF > 1%
tempdf = tempdf[tempdf['AF'].fillna(0) < 0.01]#Filter out SNPS with db VAF of > 1%
tempdf = tempdf[tempdf['EUR AF'].fillna(0) < 0.01]
tempdf = tempdf[tempdf['Global AF'].fillna(0) < 0.01]
tempdf = tempdf[tempdf['Non-Fin Eur AF'].fillna(0) < 0.01]
filtereddf.loc[len(filtereddf)] = ['>1% AF in SNP DB', len(tempdf)]


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
    

