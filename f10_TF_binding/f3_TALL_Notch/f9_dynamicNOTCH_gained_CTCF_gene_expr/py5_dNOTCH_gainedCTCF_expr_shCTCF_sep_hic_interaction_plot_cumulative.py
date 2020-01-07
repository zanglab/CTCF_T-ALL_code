import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.6)
sns.set_style("whitegrid", {'axes.grid' : False})
import CTCF_TALL_modules_new
import scipy
from scipy import stats
import bisect
sns.set_style("ticks")



def cumulative_compare_save_file(a1,b1,c1,df,xticklabels,figname):
    
    a = df.loc[a1]['log2FoldChange'].values
    b = df.loc[b1]['log2FoldChange'].values
    c = df.loc[c1]['log2FoldChange'].values
    
    data=[c,b]
    positions=[1,2]
    colors=['b','red']
#     data=[a,b]
#     positions=[0,1]
#     colors=['silver','red']
    plt.figure(figsize=(2,3))
    g1 = sns.kdeplot(c,cumulative=True,c='b')   
    g2 = sns.kdeplot(b,cumulative=True,c='red')    
    
    
    plt.ylabel('log2FC',fontsize=18)
#     plt.axes().set_ylim([0,8])
#     plt.axes().tick_params(axis='x',direction='out', length=4, width=1, colors='black')    
#     plt.axes().tick_params(axis='y',direction='out', length=4, width=1, colors='black')  
    sns.despine(offset=0, trim=False)  
#     plt.axes().set_xticklabels(xticklabels,rotation=30,fontsize=18,ha='right')
    plt.savefig(figname,bbox_inches='tight',pad_inches=.1,transparent=True,dpi=600)
    plt.close()
    

def read_genes(flag,flag2,indir):
    outfile = indir+os.sep+'{}_{}_allgenes.txt'.format(flag2,flag)
    with open(outfile,'r') as outf:
        genes = [ii.strip() for ii in outf.readlines()]
    
    outfile = indir+os.sep+'{}_{}_shCTCF_dngenes.txt'.format(flag2,flag)
    with open(outfile,'r') as outf:
        genes2 = [ii.strip() for ii in outf.readlines()]
        
    return genes, genes2


def main(outdir):

    indir = 'f3_dNOTCH_gainedCTCF_expr_in_shCTCF'
    outdir = 'f5_dNOTCH_gainedCTCF_expr_in_shCTCF_figs_cumulative'
    os.makedirs(outdir,exist_ok=True)

    deseq2_file = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/18_MYC-ChIP_shCTCF-RNA/shCTCF_RNA/f0_processing/salmon_Deseq2/salmon_Deseq2_pca/CUTLL1/f4_deseq_out_shrink/treated_CUTLL1_shCTCF_vs_ctrl_CUTLL1_PIG.csv"
    with open(deseq2_file) as inf:
        df = pd.read_csv(inf,index_col=0)
    
    flags = ['all_genes','dNotch_hic_increased_gained_ctcf_intra_domain','dNotch_hic_NOT_increased_gained_ctcf_intra_domain']
    xticklabels = ['all','dNOTCH w/ gained CTCF\n(increased Hi-C contact)']
    xticklabels = ['Decreased interaction','Increased interaction']
    
    flag2 = 'cutll1_upgenes'
    a1,a2 = read_genes(flags[0],flag2,indir)#;print(len(a1),len(a2))
    b1,b2 = read_genes(flags[1],flag2,indir)
    c1,c2 = read_genes(flags[2],flag2,indir)
    
    figname = outdir+os.sep+'dNOTCH_gained_CTCF_hic_increased_shCTCF_expr.pdf'
#     box_compare_save_file(a1,b1,c1,df,xticklabels,figname)
    cumulative_compare_save_file(a1,b1,c1,df,xticklabels,figname)
    





if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
#     parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.outdir)
