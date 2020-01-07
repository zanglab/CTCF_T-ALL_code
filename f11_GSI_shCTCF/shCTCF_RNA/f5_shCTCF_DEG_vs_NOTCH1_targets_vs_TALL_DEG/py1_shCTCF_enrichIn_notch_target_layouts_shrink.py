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
sns.set(font_scale=1.4)
sns.set_style("whitegrid", {'axes.grid' : False})
from scipy import stats
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
from matplotlib_venn import venn3,venn2
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"


def plot_scatter_ma(df,color):
    g = plt.scatter(np.log2(df['baseMean']+0.01),df['log2FoldChange'],c=color,s=23)
    return g


def plot_scatter_volcano(df,color):
    g = plt.scatter(df['log2FoldChange'],-1*np.log10(df['padj']),c=color,s=23)
    return g



outdir='f1_shCTCF_deg_vs_notch1_targets_new_figs_shrink'
os.makedirs(outdir,exist_ok=True)

fcthre,pthre=0.58,0.01
deseq2_CUTLL1_shCTCF = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/18_MYC-ChIP_shCTCF-RNA/shCTCF_RNA/f0_processing/salmon_Deseq2/salmon_Deseq2_pca/CUTLL1/f3_deseq_out/treated_CUTLL1_shCTCF_vs_ctrl_CUTLL1_PIG.csv"
deseq2_CUTLL1_shCTCF = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/18_MYC-ChIP_shCTCF-RNA/shCTCF_RNA/f0_processing/salmon_Deseq2/salmon_Deseq2_pca/CUTLL1/f4_deseq_out_shrink/treated_CUTLL1_shCTCF_vs_ctrl_CUTLL1_PIG.csv"
with open(deseq2_CUTLL1_shCTCF) as inf:
    df = pd.read_csv(inf,index_col=0)
upgenes = df[(df['log2FoldChange']>fcthre)&(df['pvalue']<pthre)].index.values
dngenes = df[(df['log2FoldChange']<-1*fcthre)&(df['pvalue']<pthre)].index.values

notch_targets_genelist = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/18_MYC-ChIP_shCTCF-RNA/shCTCF_RNA/f5_shCTCF_DEG_vs_NOTCH1_targets_vs_TALL_DEG/data/notch1_target_dynamic_RP_and_GSIwo_up.txt'
notch_targets = [i.strip() for i in open(notch_targets_genelist).readlines()]

ctrl_cols = ['CUTLL1_PIG_rep1', 'CUTLL1_PIG_rep2', 'CUTLL1_PIG_rep3']
treat_cols = ['CUTLL1_shCTCF_rep1', 'CUTLL1_shCTCF_rep2', 'CUTLL1_shCTCF_rep3']
avg_vals = (df[ctrl_cols].mean(axis=1)+df[treat_cols].mean(axis=1))*0.5
plt.figure(figsize=(3.,3.))
a = plot_scatter_ma(df,'silver')
b = plot_scatter_ma(df.loc[notch_targets],color = 'red')
# a = plot_scatter_volcano(df,'silver')
# b = plot_scatter_volcano(df.loc[notch_targets],'r')


plt.axhline(y=0,c='k',ls='--')
plt.legend([a,b],['All genes','NOTCH1 targets'],bbox_to_anchor=[-.05,1.25],loc='upper left',markerscale=1.6,borderaxespad=0,labelspacing=.1,handletextpad=0.1,frameon=False,fontsize=15)
# sns.despine(offset=0, trim=False)
plt.ylim([-3,3])
plt.ylabel('log$_2$(fold change)\n shCTCF over control',fontsize=18)
plt.xlabel('log$_2$(baseMean)',fontsize=18)
figname = outdir+os.sep+'shCTCF_DEG_down_vs_notch_target.pdf'
plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
plt.close()



