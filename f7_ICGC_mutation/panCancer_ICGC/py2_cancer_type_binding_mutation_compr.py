import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
#sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
#sns.despine(offset=0, trim=True)

import CTCF_TALL_modules_new


def main():

    indir = 'f1_mutation_around_CTCF'
    outdir = 'f2_mutation_on_cancer_specific_CTCF'
    os.makedirs(outdir,exist_ok=True)

    union_CTCF = CTCF_TALL_modules_new.return_occupancy_filtered()
    cancertype_mutation_matchness = {'BRCA':'BRCA','CRC':'COAD','LUAD':'LUAD','PRAD':'PRAD','AML':'AML','PRAD_TissueAdded':'PRAD'}
    cancertypes=['BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']
    for celltype in cancertypes:
        for expand in [9,200]:
            print(celltype,expand)
            mutation_file = indir+os.sep+'{}_mutation_on_CTCF_expand{}.csv'.format(cancertype_mutation_matchness[celltype],expand)
            mutation_event_rate_file = indir+os.sep+'{}_mutation_on_CTCF_expand{}.csv.mutation_event_rate.csv'.format(cancertype_mutation_matchness[celltype],expand)
            with open(mutation_file) as inf,open(mutation_event_rate_file) as inf2:
                df = pd.read_csv(inf,sep='\t',index_col=1)
                df2 = pd.read_csv(inf2,sep='\t',header=None,index_col=0)
        
            gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding(celltype)    
            gained_mutation = df.loc[gained_df.index]
            gained_mutation.insert(1,'id',gained_mutation.index)
            gained_mutation.to_csv(outdir+os.sep+'{}_mutation_on_{}_{}_expand{}.csv'.format(cancertype_mutation_matchness[celltype],celltype,'gained',expand),sep='\t',index=False)
            # save the mutation info (ref to alt) file, keep same format as the union one  (id in the 2nd column)
            gained_mutation_rate = df2.loc[gained_df.index]
            gained_mutation_rate.to_csv(outdir+os.sep+'{}_mutation_on_{}_{}_expand{}.mutation_rate.csv'.format(cancertype_mutation_matchness[celltype],celltype,'gained',expand),sep='\t',index=True,header=None)

            lost_mutation = df.loc[lost_df.index]
            lost_mutation.insert(1,'id',lost_mutation.index)
            lost_mutation.to_csv(outdir+os.sep+'{}_mutation_on_{}_{}_expand{}.csv'.format(cancertype_mutation_matchness[celltype],celltype,'lost',expand),sep='\t',index=False)
            lost_mutation_rate = df2.loc[lost_df.index]
            lost_mutation_rate.to_csv(outdir+os.sep+'{}_mutation_on_{}_{}_expand{}.mutation_rate.csv'.format(cancertype_mutation_matchness[celltype],celltype,'lost',expand),sep='\t',index=True,header=None)
            del df
            del df2
#             const_mutation = df.loc[union_CTCF.index]
#             const_mutation.insert(1,'id',const_mutation.index)
#             const_mutation.to_csv(outdir+os.sep+'{}_mutation_on_TALL_{}_expand{}.csv'.format(celltype,'const',expand),sep='\t',index=False)
            
        

            
            
            
            


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
  
    main()
