import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
import association_with_genes
import association_with_regions
import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')

exon_df = pd.read_csv('union_exons_overlapped.bed',sep='\t',index_col=3,header=None)
intron_df = pd.read_csv('union_introns_overlapped.bed',sep='\t',index_col=3,header=None)
promoter_df = pd.read_csv('union_promoter_overlapped.bed',sep='\t',index_col=3,header=None)

annotation_df = pd.read_csv('union_binding_mid_position.bed',sep='\t',index_col=3,header=None)
annotation_df.columns = ['chr','start','end']
annotation_df['id']=annotation_df.index


def main(infile):
    
#     promoter_index = tss_df.loc[tss_df['nearestTSS']<2500].index
#     print(len(promoter_index))

    annotation_df.loc[intron_df.index,'annotation']='Intron'
    annotation_df.loc[exon_df.index,'annotation']='Exon'
    annotation_df.loc[promoter_df.index,'annotation']='Promoter'
    annotation_df['annotation'] = annotation_df['annotation'].fillna('Distal')

    #print(annotation_df);exit()
    annotation_df.to_csv('union_binding_with_gene_annotation.csv',index=False)








if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
    parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile)
