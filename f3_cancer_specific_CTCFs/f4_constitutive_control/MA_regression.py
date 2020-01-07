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

log_para = 2

def ma_scatter_for_two(df,df2,columns,figname,color='lightgrey',yhline=False,xylim=False):
    
    fig = plt.figure(figsize=(6,6))
    #df = df.dropna()
    x,y=df[columns[0]],df[columns[1]]
    x2,y2=df2[columns[0]],df2[columns[1]]
    g1 = plt.scatter(x,y,s=15,alpha=1,c=color)
    g2 = plt.scatter(x2,y2,s=15,alpha=1,c='red')
    plt.axhline(y=0,color='gray',linestyle='--',linewidth=.7)
    if yhline:
        plt.axvline(x=0,color='gray',linestyle='--',linewidth=.7)
    plt.axis('equal')
    #plt.xlabel(xlabel,fontsize=22)
    #plt.ylabel(ylabel,fontsize=22)
    if xylim:
        x1,x2,y1,y2=xylim
        plt.xlim([x1,x2])
        plt.ylim([y1,y2])
    #plt.legend([g1,g3],['all','gained'],markerscale=3)
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1)
    plt.close()





def ma_scatter(df,columns,figname,color='lightgrey',yhline=False,xylim=False):
    
    fig = plt.figure(figsize=(6,6))
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
    plt.xlabel('avg log{}'.format(log_para),fontsize=22)
    plt.ylabel('{}'.format(columns[1]),fontsize=22)
    if xylim:
        x1,x2,y1,y2=xylim
        plt.xlim([x1,x2])
        plt.ylim([y1,y2])
    #plt.legend([g1,g3],['all','gained'],markerscale=3)
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1)
    plt.close()



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

    t,c =columns[0],columns[1]
    df = df_ori[columns]
    # log_para is set in the first row
    if log_para == 2:
        df = np.log2(df+0.01)
    elif log_para == 10:
        df = np.log10(df+0.01)
    
    df['log_avg']=0.5*(df[c]+df[t])
    df['logFC']=df[t]-df[c] 
    
    ## M=kA+C
    a,b,xmean,ymean = linear_regression(df['log_avg'].values,df['logFC'].values)
    print(a,b)
    df['reg_logFC'] = a*df['log_avg']+b
    df['adj_logFC'] = df['logFC']-df['reg_logFC']
    df[columns] = df_ori[columns]
    return df


def ma_fit_add_pvalue(df,pthre=0.05,hist_plot=False,figname=None):
    ## add pvalue and plot hist fig
    miu = np.mean(df['adj_logFC'])
    std = np.std(df['adj_logFC'])
    print('u:\t',miu,'\nstd:\t',std)
    df['pvalue'] = scipy.stats.norm(miu,std).cdf(df['adj_logFC'])
    df['pvalue'] = pd.concat([2*df['pvalue'],(1-df['pvalue'])*2],axis=1).min(axis=1)
    # plot hist graph
    left = scipy.stats.norm(miu,std).ppf(pthre*0.5)
    right = scipy.stats.norm(miu,std).ppf(1-pthre*0.5)#;print(basename,left,right)
    if hist_plot:
        plt.figure(figsize=(4,3))
        plt.hist(df['adj_logFC'],bins=80)
        plt.axes().axvline(x=left,color='red',linewidth=1)
        plt.axes().axvline(x=right,color='red',linewidth=1)
        #plt.title(basename)
        #plt.xlim([-.6,.6])
        #figname=outdir+os.sep+'hist_{}_adj_logFC.png'.format(basename)
        plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1)
        plt.close()
       
    return df

 

def main():

    outdir = 'f3_dynamic_MA_figs_csv'
    os.makedirs(outdir,exist_ok=True) 
    
    co_occurrence_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f6_Jurkat_GSI/f2_binding_signal_change_fc/f3_HC_TALL_GSI_overlap_summary/f1_HC_T-ALL_GSI_binding_occurrence/HC_TALL_GSI_binding_occurrence.csv'
    with open(co_occurrence_file) as co_occurrence_inf:
        co_occurrence_df = pd.read_csv(co_occurrence_inf,index_col = 0)    
    co_occurrence_df = co_occurrence_df.loc[CTCF_TALL_modules.return_ctcf_filtered_ids()]
    df = co_occurrence_df[(co_occurrence_df['JURKAT_DMSO']!=0)|(co_occurrence_df['JURKAT_GSI_3d']!=0)|(co_occurrence_df['JURKAT_GSI_3d_w4hr']!=0)]
    
    #### add adjusted logFC and pvalue
    basename = 'w4hr_over_gsi'
    columns = ['jurkat_gsi_3d_RPKM','jurkat_gsi_3d_w4h_RPKM']
    ma_df = ma_fit(df,columns,basename,outdir)   
    df['{}_log_avg'.format(basename)] = ma_df['log_avg'].round(3)
    df['{}_adj_logFC'.format(basename)] = ma_df['adj_logFC'].round(3)
    df['{}_adj_logFC_pvalue'.format(basename)] = ma_df['pvalue']
    
    
    #### add adjusted log FC values 
    basename = 'gsi_over_dmso'   
    columns = ['jurkat_dmso_RPKM','jurkat_gsi_3d_RPKM']
    ma_df = ma_fit(df,columns,basename,outdir)    
    df['{}_log_avg'.format(basename)] = ma_df['log_avg'].round(3)
    df['{}_adj_logFC'.format(basename)] = ma_df['adj_logFC'].round(3)
    df['{}_adj_logFC_pvalue'.format(basename)] = ma_df['pvalue']
    df.to_csv(outdir+os.sep+'Jurkat_GSI_union_bindings_MA_adj_logFC.csv')
    
   
    # MA plot
    xylim=-0.3,1.35,-.75,.75
    basename = 'w4hr_over_gsi'
    xlabel='avg logFC'
    ylabel='logFC {}'.format(basename)
    columns = ['{}_log_avg'.format(basename),'{}_adj_logFC'.format(basename)]
    figname=outdir+os.sep+'ma_{}_scatter.png'.format(basename)
    scatter_plot_for_diffGenes(df,columns,figname,xlabel,ylabel,xylim=xylim)

    # MA plot
    basename = 'gsi_over_dmso'
    columns = ['{}_log_avg'.format(basename),'{}_adj_logFC'.format(basename)]
    ylabel='logFC {}'.format(basename)
    figname=outdir+os.sep+'ma_{}_scatter.png'.format(basename)
    scatter_plot_for_diffGenes(df,columns,figname,xlabel,ylabel,xylim=xylim)

    
    # scatter plot compr log FC
    columns = ['gsi_over_dmso_adj_logFC','w4hr_over_gsi_adj_logFC']
    xlabel=columns[0]
    ylabel=columns[1]
    xylim=-.75,.75,-.75,.75
    figname=outdir+os.sep+'logFC_compr_ma.png'
    scatter_plot_for_diffGenes(df,columns,figname,xlabel,ylabel,yhline=True,xylim=xylim)




 
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

