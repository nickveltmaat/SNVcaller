import pandas as pd
import numpy as np
import sys
import io
import os
import re


sample = sys.argv[1]
data = sys.argv[2]
path = './output/'

#Load annotated variants & Filter df
annotated = pd.read_excel(path + sample + '/' + data, sheet_name='Variant', skiprows=1, engine = 'openpyxl')
if data == 'annotated_SNVs_PoN_blacklisted.xlsx':
    dffilter = pd.read_csv(path + sample + '/' +'filters_PoN.txt', sep = '\t')
elif data == 'annotated_SNVs.xlsx':
    dffilter = pd.read_csv(path + sample + '/' +'filters.txt', sep = '\t')
else: 
    print('wrong input for post filtering...')

### START FILTERING
annotated['vaf'] = annotated['Samples'].str.split('AF:', 1).str[1].str.split('_DP', 1).str[0].astype(float)
annotated['dp'] = annotated['Samples'].str.split('DP:', 1).str[1].str.split(';', 1).str[0].astype(int)
annotated['mut_dp'] = (annotated['vaf'] * annotated['dp']).round().astype(int)
annotated = annotated[annotated['mut_dp'] >= 2] #Filter out variants with 1 mutant read
dffilter.loc[len(dffilter)] = ['<2_mut_reads', len(annotated)]
annotated = annotated[annotated['mut_dp'] >= 3] #Filter out variants with 2 mutant reads
dffilter.loc[len(dffilter)] = ['<3_mut_reads', len(annotated)]

#drop created new colums
annotated = annotated.drop(['vaf', 'dp', 'mut_dp'], axis = 1)

#annotated['AF'] = annotated['AF'].fillna(0)
#annotated = annotated[annotated['AF'] < 0.01] #Filter 1000G
#annotated['AF'] = annotated['AF'].replace(0, 'NaN', inplace=True)
#dffilter.loc[len(dffilter)] = ['>1%_1000G', len(annotated)]

#annotated['EUR AF'] = annotated['EUR AF'].fillna(0)
#annotated = annotated[annotated['EUR AF'] < 0.01] #Filter 1000G
#annotated['EUR AF'] = annotated['EUR AF'].replace(0, 'NaN', inplace=True)
#dffilter.loc[len(dffilter)] = ['>1%_1000G_EUR', len(annotated)]

#annotated['Global AF'] = annotated['Global AF'].fillna(0)
#annotated = annotated[annotated['Global AF'] < 0.01] #Filter GnoMAD
#annotated['Global AF'] = annotated['Global AF'].replace(0, 'NaN', inplace=True)
#dffilter.loc[len(dffilter)] = ['>1%_GnomAD', len(annotated)]

#annotated['Non-Fin Eur AF'] = annotated['Non-Fin Eur AF'].fillna(0)
#annotated = annotated[annotated['Non-Fin Eur AF'] < 0.01] #Filter GnoMAD
#annotated['Non-Fin Eur AF'] = annotated['Non-Fin Eur AF'].replace(0, 'NaN', inplace=True)
#dffilter.loc[len(dffilter)] = ['>1%_GnomAD_EUR', len(annotated)]

annotated = annotated[annotated['Sequence Ontology'] != 'synonymous_variant'] #Filter synonymous
dffilter.loc[len(dffilter)] = ['Synonymous', len(annotated)]

annotated = annotated[annotated['Sequence Ontology'] != '3_prime_UTR_variant'] #Filter 3' UTR variants
dffilter.loc[len(dffilter)] = ["3' UTR", len(annotated)]

annotated = annotated[annotated['Sequence Ontology'] != '5_prime_UTR_variant'] #Filter 5' UTR variants
dffilter.loc[len(dffilter)] = ["5' UTR", len(annotated)]

annotated = annotated[annotated['Sequence Ontology'] != 'intron_variant'] #Filter Intronic variants
dffilter.loc[len(dffilter)] = ['Intronic', len(annotated)]

annotated = annotated[annotated['Sequence Ontology'] != '2kb_upstream_variant'] #Filter 2kb Upstream variants
dffilter.loc[len(dffilter)] = ['2kb upstream', len(annotated)]

annotated = annotated[annotated['Sequence Ontology'] != '2kb_downstream_variant'] #Filter 2kb Downstream variants
dffilter.loc[len(dffilter)] = ['2kb downstream', len(annotated)]


#Save Filtered variants & Filtering info
if data == 'annotated_SNVs_PoN_blacklisted.xlsx':
    dffilter.to_csv(path + sample + '/filters_PoN.txt', sep='\t', index = False, header= True)
    annotated.to_excel(path + sample + '/filtered_variants_PoN.xlsx', index=False)
    
elif data == 'annotated_SNVs.xlsx':
    dffilter.to_csv(path + sample + '/filters.txt', sep='\t', index = False, header= True)
    annotated.to_excel(path + sample + '/filtered_variants.xlsx', index=False)
else: 
    print('wrong input for post filtering...')
