'''
select within domain interactions between jurkat and cd4 then cal logFC
'''

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
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
import association_with_genes
import association_with_regions
import re,bisect
import CTCF_TALL_modules
import scipy
from scipy import stats
import MA_regression

def linear_regression(x,y):
    xmean = np.mean(x)
    ymean = np.mean(y)
    assert len(x)==len(y)
    sum1,sum2=0,0
    for i in np.arange(len(x)):
        sum1 += (x[i]-xmean)*(y[i]-ymean)
        sum2 += np.power((x[i]-xmean),2)
    a = sum1/sum2
    b = ymean-a*xmean
    return a,b,xmean,ymean



def ma_fit(df_ori,columns):

    t,c = columns[0],columns[1]
    df = df_ori[columns]
    df = np.log2(df+0.001)
    
    df['log_avg']=0.5*(df[c]+df[t])
    df['logFC']=df[t]-df[c] 
    
    ## M=kA+C
    a,b,xmean,ymean = linear_regression(df['log_avg'].values,df['logFC'].values)
    print(a,b)
    df['reg_logFC'] = a*df['log_avg']+b
    df['adj_logFC'] = df['logFC']-df['reg_logFC']
    df[columns] = df_ori[columns]
    return df


def gained_lost_scatter(x,y,figname,cancertype,xlabel,ylabel,hline):
    
    fig = plt.figure(figsize=(6,6))
    g1 = plt.scatter(x,y,c='lightgrey',s=15,alpha=1,cmap = plt.cm.GnBu_r)
    
    gained_df,lost_df,const_df = CTCF_TALL_modules.return_specific_binding_df(cancertype)
    const_df = pd.read_csv('/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f2_bindings_signals_fdr/f3_proper_ctrl_constitutive_selection/const_NONoverlapped_panCancer_fdr_LT_001_TALL_log2fcLT1.bed',sep='\t',index_col=3,header=None)
    constitutive = CTCF_TALL_modules.return_constitutive_20k_df()
    #print(x[gained_df.index]);exit()
    g4 = plt.scatter(x[constitutive.index],y[constitutive.index],c='lightgreen',s=15,alpha=.5,cmap = plt.cm.GnBu_r)
    g5 = plt.scatter(x[const_df.index],y[const_df.index],c='green',s=15,alpha=.5,cmap = plt.cm.GnBu_r)
    g2 = plt.scatter(x[gained_df.index],y[gained_df.index],c='r',s=15,alpha=.5,cmap = plt.cm.GnBu_r)
    g3 = plt.scatter(x[lost_df.index],y[lost_df.index],c='b',s=15,alpha=.5,cmap = plt.cm.GnBu_r)
    plt.axhline(y=hline,color='gray',linestyle='--',linewidth=.7)
    plt.title(cancertype)
    plt.xlabel(xlabel,fontsize=22)
    plt.ylabel(ylabel,fontsize=22)
    #plt.xlim([-5,5])
    #plt.ylim([y1,y2])
    plt.legend([g1,g2,g3,g4,g5],['all','gained','lost','constitutive','const-ctrl'],markerscale=2,fontsize=12)
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1)
    plt.close()



def cancertype_maplot(cancertype,compr_type):

    outdir = 'f5_py3_ma_plot'
    os.makedirs(outdir,exist_ok=True) 

    infile = 'f1_cancertype_ttest_fdr_logFC/{}_signals_ttest_fdr.csv'.format(cancertype)
    with open(infile) as inf:
        df = pd.read_csv(inf,index_col =0,sep='\t',low_memory=False)    
    #print(df.columns);exit()
    
    # scatter plot of logFC vs. adjp from ttest on signals
    df['cancer_ctrl_fc'] = np.log2(df['cancer_mean']+0.001)-np.log2(df['{}_mean'.format(compr_type)]+0.001)
    figname = outdir+os.sep+'valcano_cancer_vs_{}_{}.png'.format(compr_type,cancertype)
    x,y = df['cancer_ctrl_fc'],-np.log10(df['cancer_vs_{}_adjp'.format('normal' if compr_type=='normal' else 'other')])
    #density_scatter(x,y,figname,cancertype,color='lightgrey',yhline=False,xylim=False)
    xlabel,ylabel = 'log2 FC','-log10 adjp'
    gained_lost_scatter(x,y,figname,cancertype,xlabel,ylabel,hline=-np.log10(0.05))
    
    
    # plot MA plot of binding signals between cancer vs. ctrl
    figname = outdir+os.sep+'MA_cancer_vs_{}_{}.png'.format(compr_type,cancertype)
    compr_columns = ['cancer_mean','{}_mean'.format(compr_type)]
    df = ma_fit(df,compr_columns)
    #gained_df,lost_df,const_df = CTCF_TALL_modules.return_specific_binding_df(cancertype)
    #print(df.loc[gained_df.index]);exit()
    #a = df.loc[df_fit[df_fit['logFC']>3.5].index]
    #a = a[a['cancer_mean']>10]
    #print(a);exit()
    x,y = df['log_avg'],df['logFC']
    xlabel,ylabel =  'log2_avg','log2FC'
    gained_lost_scatter(x,y,figname,cancertype,xlabel,ylabel,hline=0)
    #exit()
    


def main():

    for cancertype in ['T-ALL','Breast_cancer','Colon_cancer','Lung_cancer']:
        for compr_type in ['normal','ctrl']:
            cancertype_maplot(cancertype,compr_type)





 
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

