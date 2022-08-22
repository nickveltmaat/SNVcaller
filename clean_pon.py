from functionsmodule import *

inputfile = sys.argv[1]
outputfile = sys.argv[2]

#Load df, add zeroes to tools called and remove Chr
ufblacklist = pd.read_excel(inputfile, engine='openpyxl', sheet_name = 'Variant', skiprows = 1)

ufblacklist['Samples'] = ufblacklist['Samples'].apply(lambda x: str(x).zfill(4))
ufblacklist['Chrom'] = ufblacklist['Chrom'].str.replace('chr', '')

#start filtering here
blacklist = ufblacklist.copy()
#blacklist['Clinical Significance'] = blacklist['Clinical Significance'].fillna('NaN')
#blacklist = blacklist[~blacklist['Clinical Significance'].str.contains("atho")] #Remove P/LP vars

#Remove Synonymous, UTR, up-/downstream and intronic variants
blacklist = blacklist[~blacklist['Sequence Ontology'].str.contains("synonymous")]
blacklist = blacklist[~blacklist['Sequence Ontology'].str.contains("UTR")]
blacklist = blacklist[~blacklist['Sequence Ontology'].str.contains("kb_")]
blacklist = blacklist[~blacklist['Sequence Ontology'].str.contains("intron")]

#Final blacklist: Select columns
df = blacklist[["Chrom", "Pos", "Reference allele", "Alternate allele", "Samples"]].copy()

df.to_csv(outputfile, header=False, index=False, sep='\t')
