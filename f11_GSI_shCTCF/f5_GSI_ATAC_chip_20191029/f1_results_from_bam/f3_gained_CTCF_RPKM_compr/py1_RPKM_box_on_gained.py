import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=17
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
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"





def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_rel(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),100)*1.01 ,1.02, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),0)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='*'
#     p_label='{:.2e}'.format(p)
#     if p_label[-2]=='0':
#         p_label = p_label[:-2]+p_label[-1]
    if p<1:
        if p<0.001:
            p_label = '**'
        if compr_pos[2] == 't':
            plt.plot([x1*1.03, x1*1.03, x2*0.95, x2*0.95], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*h, p_label, ha='center', va='bottom', color=col,fontsize=18)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.95, x2*0.95], [y2, y2*.91, y2*.91, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*.95, p_label, ha='center', va='top', color=col,fontsize=18)



def box_plot(df,gained_df,columns,outdir,suffix):
    
    df = df.loc[gained_df.index].dropna() ;print(df.shape)
#     df = np.log2(df+1)
    fig = plt.figure(figsize=(2.6,2.6))
    positions = [0,1,2]
    colors = ['grey','grey','grey']
    box_vals = [df[columns[0]],df[columns[1]],df[columns[2]]]

    g = plt.boxplot(box_vals,positions=positions,widths = .6,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='silver'),showfliers=False)    
    for patch, color in zip(g['boxes'], colors):
        patch.set_facecolor(color)

    scatter_X = []
    for position_id in np.arange(len(positions)):
        scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
        plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=20,zorder=0,alpha=0.99)

    positions=[0,1,2]
    for compr_pos in [[0,1,'t'],[1,2,'t']]:
        mark_pvalue(compr_pos,positions,box_vals)

#     sns.despine(offset=0, trim=False)
    plt.axes().set_xticks(positions)
    plt.axes().set_xticklabels(['DMSO','GSI-3d','GSI-washout'],rotation=30,ha='right') 
#     plt.axhline(y=0,color='grey',linestyle='--',linewidth=.7)
#     plt.legend([g["boxes"][0],g["boxes"][1]],['Control','Gained'],bbox_to_anchor=[1,1])
    plt.ylabel('ATAC-seq RPKM')
    plt.ylim([-1,19])
    #plt.xlabel(columns[0])
    #plt.title('{} dynamic change'.format(hm_name))
    figname = outdir+os.sep+'{}_box.pdf'.format(suffix)
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,transparent=True)
    plt.close()


  

def main():

    outdir = 'f1_box_figs'
    os.makedirs(outdir,exist_ok=True)
    
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    const_df = CTCF_TALL_modules_new.return_constitutive_df()
    gained_df = gained_df[gained_df['motif_strand']!='N']
    
    indir='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f11_GSI_shCTCF/f5_GSI_ATAC_chip_20191029/f1_results_from_bam/f1_binding_signal_csv/binding_changes'
    suffixs = ['e200','e150','binding_site']
#     suffixs = ['peak_overlapped_union_e200','peak_overlapped_union_e500','peak_overlapped_union_binding_site']
    treats = ['dmso','gsi_3d','gsi_wo']

    for suffix in suffixs:
        df_file=indir+os.sep+'ATAC_RPKM_changes_on_Union_{}.csv'.format(suffix)
        df = pd.read_csv(df_file,index_col=0,sep=',')#;print(df)
        
        box_plot(df,gained_df,treats,outdir,suffix)
        
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
