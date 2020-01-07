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
# import CTCF_TALL_modules
import scipy
from scipy import stats



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



def ma_fit(df,columns,basename,outdir,pthre=0.05):

    c,t=columns[0],columns[1]
    df = df[columns]
    df = np.log2(df+0.01)
    
    df['log_avg']=0.5*(df[c]+df[t])
    df['logFC']=df[t]-df[c] 
    
    ## M=kA+C
    a,b,xmean,ymean = linear_regression(df['log_avg'].values,df['logFC'].values)
    print(a,b)
    df['reg_logFC'] = a*df['log_avg']+b
    df['adj_logFC'] = df['logFC']-df['reg_logFC']

    ## add pvalue and plot hist fig
    miu = np.mean(df['adj_logFC'])
    std = np.std(df['adj_logFC'])
    print('u:\t',miu,'\nstd:\t',std)
    df['pvalue'] = scipy.stats.norm(miu,std).cdf(df['adj_logFC'])
    df['pvalue'] = pd.concat([2*df['pvalue'],(1-df['pvalue'])*2],axis=1).min(axis=1)
    # plot hist graph
    left = scipy.stats.norm(miu,std).ppf(pthre*0.5)
    right = scipy.stats.norm(miu,std).ppf(1-pthre*0.5);print(basename,left,right)
    plt.figure(figsize=(4,3))
    plt.hist(df['adj_logFC'],bins=80)
    plt.axes().axvline(x=left,color='red',linewidth=1)
    plt.axes().axvline(x=right,color='red',linewidth=1)
    plt.title(basename)
    #plt.xlim([-.6,.6])
    figname=outdir+os.sep+'hist_{}_adj_logFC.png'.format(basename)
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1)
    plt.close()
       
    return df

 

def main():

    outdir = 'f2_dynamic_MA_csv_log2'
    os.makedirs(outdir,exist_ok=True) 
    
    df = pd.DataFrame()
    file_pardir = 'f1_signal_csv'    
    for gsi_index in ['jurkat_dmso','jurkat_gsi_3d','jurkat_gsi_3d_w4h']:
        gsi_index_file = file_pardir+os.sep+'{}.csv'.format(gsi_index)
        with open(gsi_index_file) as gsi_index_inf:
            gsi_index_df = pd.read_csv(gsi_index_inf,sep='\t',index_col=0,header=None)
        gsi_index_df.columns = ['{}_RPKM'.format(gsi_index)]
        gsi_index_df = gsi_index_df.round(2)#;print(gsi_index_df)
        df = pd.concat([df,gsi_index_df[['{}_RPKM'.format(gsi_index)]]],axis=1)
#     print(df);exit()
    
    # keep only bindings show at least one peak in GSI ChIP
    gsi_peak_overlapped_infile = 'f0_GSI_peak_overlapped_union_binding/union_binding_GSI_peak_overlapped.bed'
    with open(gsi_peak_overlapped_infile) as gsi_peak_overlapped_inf:
        df_peak_overlapped = pd.read_csv(gsi_peak_overlapped_inf,sep='\t',header=None,index_col=3)    
    
    df = df.loc[df_peak_overlapped.index]
   
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

