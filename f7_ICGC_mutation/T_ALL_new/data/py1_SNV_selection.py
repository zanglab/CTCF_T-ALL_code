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

from collections import Counter


def main():

    indir = 'ctcf_sites_output'
    outdir = 'ctcf_sites_output_SNPs'
    os.makedirs(outdir,exist_ok=True)
    
    suffix='.tbl'
    infiles = glob.glob(indir+os.sep+'*{}'.format(suffix))
    mutation_counter = set()
    for infile in infiles:
        basename = os.path.basename(infile).split(suffix)[0]
        with open(infile) as inf:
            df = pd.read_csv(inf,sep='\t',index_col=None,header=None)
        df.columns = ['chr','pos','ref','alt','genotype']
        df['ref_len']=[len(i) for i in df['ref']]
        df['alt_len']=[len(i) for i in df['alt']]
        df = df[(df['ref_len']==1)&(df['alt_len']==1)]
        df = df[['chr','pos','ref','alt','genotype']]
        mutation_counter = Counter(df['ref']+'>'+df['alt'])
        outfile1 = '{}/{}{}'.format(outdir,basename,suffix)
        outfile2 = '{}/{}{}.counter'.format(outdir,basename,suffix)
        df.to_csv(outfile1,sep='\t',index=False)
        with open(outfile2,'w') as outf2:
            for key in sorted(mutation_counter.keys(),key=lambda i: mutation_counter[i],reverse=True):
                outf2.write('{}\t{}\n'.format(key,mutation_counter[key]))
        #print(df);exit()



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
