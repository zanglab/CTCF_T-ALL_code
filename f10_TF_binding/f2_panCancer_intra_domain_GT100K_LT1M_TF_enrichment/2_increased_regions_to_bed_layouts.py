import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
# import association_with_genes
# import association_with_regions
# import re,bisect
# import CTCF_TALL_modules
import scipy
from scipy import stats
#import MA_regression


def write_out_bed(df,outfile):

    with open(outfile,'w') as outf:
        for i in df.index:
            a,b,c = i.split('_')
            fc = df.loc[i,'adj_logFC']
            pos = int(b)+int(c)#;print(a,pos);exit()
            outf.write('{}\t{}\t{}\t{}\n'.format(a,pos-5000,pos,fc))
    outf.close()


def ma_scatter(df,columns,figname,color='lightgrey',yhline=False,xylim=False):
    
    fig = plt.figure(figsize=(3.5,3.5))
    #df = df.dropna()
    x,y=df[columns[0]],df[columns[1]]
    xy = np.vstack([x,y])
    z = scipy.stats.gaussian_kde(xy)(xy)
    #plt.scatter(x,y,c=z,s=5,cmap = plt.cm.GnBu_r)
    g1 = plt.scatter(x,y,c=z,s=15,alpha=1,cmap = plt.cm.GnBu_r)
    plt.axhline(y=0,color='gray',linestyle='--',linewidth=.7)
    if yhline:
        plt.axvline(x=0,color='gray',linestyle='--',linewidth=.7)
    plt.axis('equal')
    plt.xlabel('Average (log2 Hi-C interaction)',fontsize=22)
    plt.ylabel('log2 FC',fontsize=22)
    if xylim:
        x1,x2,y1,y2=xylim
        plt.xlim([x1,x2])
        plt.ylim([y1,y2])
    #plt.legend([g1,g3],['all','gained'],markerscale=3)
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,dpi=600,transparent=True)
    plt.close()



def retion_to_bed(df,outdir,basename):

    columns = df.columns[:2]
    value_df = df[columns]
    
    side_index = value_df.loc[(value_df.T==0).any()].index
    side_df = df.loc[side_index]
    center_index = value_df.loc[(value_df.T!=0).all()].index
    center_df = df.loc[center_index]
    
    #print(center_df);exit()
    
    side_hicor = side_df.loc[side_df['logFC']>np.log2(300)];print(basename,'\nside',side_hicor.shape)
    side_ctrl = side_df.loc[(side_df['logFC']<np.log2(60))&(side_df['logFC']>0)];print('side-ctrl',side_ctrl.shape)

    center_hicor = center_df.loc[(center_df['logFC']>np.log2(2))&(center_df['log_avg']>0.05)];print('center',center_hicor.shape)
    center_ctrl = center_df.loc[(center_df['logFC']>0)&(center_df['logFC']<0.5)&(center_df['log_avg']<0.05)];print('center-ctrl',center_ctrl.shape)
    ctrl_df = pd.concat([side_ctrl,center_ctrl])
    #print(ctrl_df)
    outname = outdir+os.sep+'{}_{}'.format(basename,'side_hicor')
    columns = ['log_avg','logFC']
    ma_scatter(df,columns,'{}.pdf'.format(outname))
    #write_out_bed(side_hicor,'{}.bed'.format(outname))
    

    

def main():

    outdir = 'f2_high_cor_high_logFC_increased_layouts'
    os.makedirs(outdir,exist_ok=True) 


    for basename in ['T-ALL_gained']:
	    logFC_file = 'f1_high_cor_high_logFC_regions/{}_logFC.bed'.format(basename)
	
	    with open(logFC_file) as logFC_inf:
		    logFC_df = pd.read_csv(logFC_inf,index_col=0)
		    logFC_df = logFC_df.dropna(how='all')
		    retion_to_bed(logFC_df,outdir,basename)


 
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

