import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=12
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
import scipy
from scipy import stats
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
import CTCF_TALL_modules_new



def plot_each_point(df,df_id,color,size=3):
    
    plot_df = df.loc[df_id.index].dropna()
    g = plt.scatter(plot_df[plot_df.columns[0]],plot_df[plot_df.columns[1]],c = color,s=size)
    return g

    
def scatter_plot(df,gained_df,lost_df,const_df,columns,outdir,suffix):

    df = df[columns]
    figname = outdir+os.sep+'{}_scatter.png'.format(suffix)
    
    fig = plt.figure(figsize=(4,4))
    plt.axhline(y=0,color='grey',linestyle='--',linewidth=.5)
    plt.axvline(x=0,color='grey',linestyle='--',linewidth=.5)
    
    a = plot_each_point(df,const_df,'lightgrey',size=5) 
    b = plot_each_point(df,gained_df,'red',size=5) 
    c = plot_each_point(df,lost_df,'blue',size=5) 

    plt.legend([a,b,c],['Control','Gained','lost'], markerscale=3,fontsize=10,loc=0,frameon=False)
    plt.xlabel(columns[0])
    plt.ylabel(columns[1])
    #plt.axis('equal')
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,transparent=True)
    plt.close()


def mark_pvalue(compr_pos,positions,box_vals):

    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),99)*1.3 ,1.05, 'k'
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    plt.plot([x1, x1, x2, x2], [y, y*h, y*h, y], lw=1, c=col)
    if p<1:
        plt.text((x1+x2)*.5, y*h, "{:.2f}".format(p), ha='center', va='bottom', color=col,fontsize=18)
    else:
        plt.text((x1+x2)*.5, y*h, "n.s.", ha='center', va='bottom', color=col,fontsize=18)



def box_plot(df,gained_df,lost_df,const_df,columns,outdir,suffix):
    
    df1 = df.loc[const_df.index].dropna() ;print(df1.shape)# const
    df2 = df.loc[gained_df.index].dropna() ;print(df2.shape)# gained
    
    figname = outdir+os.sep+'{}_box.png'.format(suffix)
    
    fig = plt.figure(figsize=(3,3))
    positions = [1,1.5,2.5,3]
    
    box_vals = [df1[columns[0]],df2[columns[0]],df1[columns[1]],df2[columns[1]]]
    g = plt.boxplot(box_vals,positions=positions,widths = .4,patch_artist=True,showfliers=False)
    
    compr_pos=[0,1]
    mark_pvalue(compr_pos,positions,box_vals)

    compr_pos=[2,3]
    mark_pvalue(compr_pos,positions,box_vals)

    colors = [ 'lightgrey','red', 'lightgrey','red']
    for patch, color in zip(g['boxes'], colors):
        patch.set_facecolor(color)
    sns.despine(offset=0, trim=False)
    plt.axes().set_xticks(positions[0:-1:2])
    plt.axes().set_xticklabels(columns,rotation=15,ha='center') 
    plt.axhline(y=0,color='grey',linestyle='--',linewidth=.7)
    plt.legend([g["boxes"][0],g["boxes"][1]],['Control','Gained'],bbox_to_anchor=[1,1])
    #plt.ylabel(columns[1])
    #plt.xlabel(columns[0])
    #plt.title('{} dynamic change'.format(hm_name))
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1)
    plt.close()


  

def main():

    outdir = 'f3_scatter_plot_figs_col_row_norm'
    os.makedirs(outdir,exist_ok=True)
    
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    const_df = CTCF_TALL_modules_new.return_constitutive_df()
    
    indir='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f11_GSI_shCTCF/f5_GSI_ATAC_chip_20191029/f1_results_from_bam/f1_binding_signal_csv/binding_changes_col_row_norm'
    suffixs = ['e200','e500','binding_site','peak_overlapped_union_e200','peak_overlapped_union_e500','peak_overlapped_union_binding_site']
#     suffixs = ['peak_overlapped_union_e200','peak_overlapped_union_e500','peak_overlapped_union_binding_site']
    treats = ['dmso','gsi_3d','gsi_wo']

    for suffix in suffixs:
        df_file=indir+os.sep+'ATAC_RPKM_changes_on_Union_{}.csv'.format(suffix)
        df = pd.read_csv(df_file,index_col=0,sep=',')#;print(df)
        
        columns = ['gsi_3d_vs_dmso','gsi_wo_vs_gsi_3d']
        scatter_plot(df,gained_df,lost_df,const_df,columns,outdir,suffix)            
        box_plot(df,gained_df,lost_df,const_df,columns,outdir,suffix)
        
#         exit()




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
