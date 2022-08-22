import pandas as pd
import numpy as np
import sys
import io
import os
import math
import glob

#Reading .vcf into pandas df:
def read_vcf(path):
    '''
    This functions reads a .vcf file and returns a pandas dataframe
    '''
    with open(path, 'r') as f:
      lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(io.StringIO(''.join(lines)), 
                       dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 
                              'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str}, 
                       sep='\t').rename(columns={'#CHROM': 'CHROM'})

#Generating RD, VAF & MDP columns from .vcf
def vcf_info_vardict(df):
    '''
    This functions loads a Vardict .vcf and extracs 'RD', 'VAF' & 'MDP' data into new columns
    '''
    df['MDP'] = df['INFO'].str.split('ADP=',1).str[1].str.split(';RFW',1).str[0].astype(int)
    df['RD'] = df['INFO'].str.split('DP=',1).str[1].str.split(';ADP',1).str[0].astype(int)
    df['AF'] = df['INFO'].str.split('AF=',1).str[1].str.split(';BIAS',1).str[0].astype(float)
    return df

def vcf_info_lofreq(df):
    '''
    This functions loads a Lofreq .vcf and extracs 'RD', 'VAF' & 'MDP' data into new columns
    '''
    df['AF'] = df['INFO'].str.split('AF=',1).str[1].str.split(';SB',1).str[0].astype(float)
    df['RD'] = df['INFO'].str.split('DP=',1).str[1].str.split(';AF',1).str[0].astype(int)
    df['MDP'] = (df['AF'] * df['RD']).round().astype(int)
    return df

def vcf_info_mutect(df):
    '''
    This functions loads a Mutect2 .vcf and extracs 'RD', 'VAF' & 'MDP' data into new columns
    '''
    df = df[df['FILTER'] ==  'PASS'].copy() #Already filter on PASS here
    df['AF'] = df[df.columns[9]].str.split(':', 5).str[2].astype(float)
    df['RD'] = df[df.columns[9]].str.split(':', 5).str[3].astype(int)
    df['MDP'] = (df['AF'] * df['RD']).round().astype(int)
    return df

def vcf_info_sinvict(df):
    '''
    This functions loads a Sinvict .vcf and extracs 'RD', 'VAF' & 'MDP' data into new columns
    '''
    df['MDP'] = df['calls_level1_sorted.sinvict'].str.split(':',2).str[2].astype(int)
    df['RD'] = df['calls_level1_sorted.sinvict'].str.split(':',2).str[0].astype(int)
    df['AF'] = df['calls_level1_sorted.sinvict'].str.split(':',2).str[1].astype(float)
    return df


def roundup(x):
    """
    This function rounds up a number to the first higher 200-fold

    Argument: Float or Int

    Returns 200-fold int
    """
    return int(math.ceil(x / 200.0)) * 200
