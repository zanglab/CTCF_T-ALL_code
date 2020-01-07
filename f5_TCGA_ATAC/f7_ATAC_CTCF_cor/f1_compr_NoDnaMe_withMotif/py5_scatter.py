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

def scatter(ctcf_df,atac_df,figname):

    fig = plt.figure(figsize=(4,4))

    x,y = ctcf_df[1].values,atac_df[1].values
    
    data, x_e,y_e = np.histogram2d(x,y,bins=20)
    z = interpn(( 0.5*(x_e[1:]+x_e[:-1]),0.5*(y_e[1:]+y_e[:-1]) ),data,np.vstack([x,y]).T, method = 'splinef2d', bounds_error=False )
    idx = z.argsort()
    x,y,z = x[idx],y[idx],z[idx]
    g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=1,marker='o')

#     plt.scatter(ctcf_df.iloc[:,:],atac_df.iloc[:,:],s=1,c='k')
    plt.xlabel('CTCF RPKM',fontsize=18)
    plt.ylabel('ATAC RPKM',fontsize=18)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()


    
def main():

    
    ctcf_file='f4_patterns/ctcf_RPKM.csv'
    atac_file='f4_patterns/atac_RPKM.csv'
    
    with open(ctcf_file) as ctcf_inf, open(atac_file) as atac_inf:
        ctcf_df = pd.read_csv(ctcf_inf,sep='\t',index_col=0,header=None)
        atac_df = pd.read_csv(atac_inf,sep='\t',index_col=0,header=None)
    ctcf_df = np.log10(ctcf_df+.1)
    atac_df = np.log10(atac_df+.1)
    figname = 'f5_figs/ctcf_cor_atac_scatter_binding_position.png'
    scatter(ctcf_df,atac_df,figname)
    
    
    ctcf_file='f4_patterns/ctcf_RPKM_e200.csv'
    atac_file='f4_patterns/atac_RPKM_e200.csv'
    
    with open(ctcf_file) as ctcf_inf, open(atac_file) as atac_inf:
        ctcf_df = pd.read_csv(ctcf_inf,sep='\t',index_col=0,header=None)
        atac_df = pd.read_csv(atac_inf,sep='\t',index_col=0,header=None)
    ctcf_df = np.log10(ctcf_df+.1)
    atac_df = np.log10(atac_df+.1)
    figname = 'f5_figs/ctcf_cor_atac_scatter_e200.png'
    scatter(ctcf_df,atac_df,figname)
    
    
            
            
            
        

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

