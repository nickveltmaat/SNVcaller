'''
This module filters Mutect2, Vardict, Lofreq or Sinvict .vcf files on:
  - 'PASS' (Mutect2)
  - Read Depth (Vardict)
  - Mutant read depth - # of reads supporting a variant (All tools)
'''
from functionsmodule import *

min_mrd = int(sys.argv[1])
stype = sys.argv[2]
tool = sys.argv[3]
output = sys.argv[4]
min_rd = int(sys.argv[5])

#Filter .vcf based on RD & MDP    
def filter_vcf(df):
    '''
    This functions loads a .vcf converted to a dataframe using read_vcf with added 'RD', 'VAF' & 'MDP' columns
    '''
    df = df[df['RD'] >= min_rd] #Filter out variants with x or less total reads
    df = df[df['MDP'] >= min_mrd].drop(['RD', 'AF', 'MDP'], axis = 1)  #Filter out variants with 2 or less mutant reads
    return df

#Writing .vcf function for Lofreq & Vardict
def write_lf_vd(description, filtereddf, outfile):
    for i in description:
        outfile.write('##'+i)
    outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t"+"\n")
    for i in filtereddf.values:
        outfile.write(i[0] +'\t'+str(i[1])+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]+'\t'+i[6]+'\t'+i[7]+'\n')
    return outfile

#Main function:    
def process_vcf(tool, **kwargs):
    snumber = kwargs.get('number', None)
    #Picking correct imput file    
    if stype == 'PoN':
        vcffile = './PoN/normal_' + snumber + '_' + tool + '.vcf'
    elif stype == 'sample':
        if tool == 'Vardict':
            vcffile = './temp/VD/VD_1.vcf'
        elif tool == 'Sinvict':
            vcffile = './temp/SV/SV_1.vcf'
        elif tool == 'Lofreq':
            vcffile = './temp/LF/LF_1.vcf'
        elif tool == 'Mutect':
            vcffile = './temp/M2/M2_1.vcf'
    
    ###Save description (every line starting with '##') in list: description
    with open(vcffile, "r") as fi:
        description = []
        for ln in fi:
            if ln.startswith("##"):
                description.append(ln[2:])
                
    ###Load the rest (not starting with ##) from vcf as pandas df
    df = read_vcf(vcffile)
    
    ###Filter dataframe & Write new .vcf file
    if tool == 'Vardict':
        filtereddf = filter_vcf(vcf_info_vardict(df))
        outfile = open(output, "w")
        write_lf_vd(description, filtereddf, outfile)
    
    elif tool == 'Lofreq':
        filtereddf = filter_vcf(vcf_info_lofreq(df))
        outfile = open(output, "w")
        write_lf_vd(description, filtereddf, outfile)
    
    elif tool == 'Mutect': 
        filtereddf = filter_vcf(vcf_info_mutect(df))
        snumber = filtereddf.columns[-1]
        outfile = open(output, "w")
        for i in description:
            outfile.write('##'+i)
        outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+snumber+"\n")
        for i in filtereddf.values:
            outfile.write(i[0] +'\t'+str(i[1])+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]+'\t'+i[6]+'\t'+i[7]+'\t'+i[8]+'\t'+i[9]+'\n')               
                       
    elif tool == 'Sinvict': 
        filtereddf = filter_vcf(vcf_info_sinvict(df))
        outfile = open(output, "w")
        for i in description:
            outfile.write('##'+i)
        outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcalls_level1_sorted.sinvict"+"\n")
        for i in filtereddf.values:
            outfile.write(i[0] +'\t'+str(i[1])+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]+'\t'+i[6]+'\t'+i[7]+'\t'+i[8]+'\t'+i[9]+'\n')
        
    outfile.close()
    print(vcffile + "  -->  " + output)
    return None

#Main calling function
if stype == 'PoN':
    snumber = sys.argv[6] #Only used for healthy/normal samples
    process_vcf(tool, number=snumber)
elif stype == 'sample':
    process_vcf(tool)
