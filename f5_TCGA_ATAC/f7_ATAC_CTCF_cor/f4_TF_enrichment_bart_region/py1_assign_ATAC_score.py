import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import gridspec
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=12
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.6)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
import CTCF_TALL_modules_new
import scipy
from scipy.interpolate import interpn

  
    
def main():        
            

    indir = '../f3_TF_enrichment_analysis/data/'
    outdir = 'data'
    os.makedirs(outdir,exist_ok=True)

    atac_file='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f5_TCGA_ATAC/f7_ATAC_CTCF_cor/f0_union_binding_RPKM_pattern/binding_csv/atac_RPKM_union_site.csv'
    atac_df = pd.read_csv(atac_file,sep='\t',index_col=0,header=None)
    
    infiles = glob.glob(indir+'/union*bed') 
    for infile in infiles:
        basename = os.path.basename(infile)
        df = pd.read_csv(infile,sep='\t',index_col=3,header=None)
        df.insert(3,'score',atac_df.loc[df.index][1])
        df.insert(3,'id',df.index)
        df.to_csv(outdir+os.sep+basename,sep='\t',header=None,index=False)

    infiles = glob.glob(indir+'/*csv') 
    for infile in infiles:
        basename = os.path.basename(infile).split('.csv')[0]
        df = pd.read_csv(infile,sep='\t',index_col=0,header=None)
        df.insert(3,'score',atac_df.loc[df.index][1])
        df.insert(3,'id',df.index)
        df.to_csv(outdir+os.sep+basename+'.bed',sep='\t',header=None,index=False)











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

