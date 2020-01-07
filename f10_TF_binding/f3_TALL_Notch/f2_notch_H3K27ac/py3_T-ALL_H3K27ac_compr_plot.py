import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
#import association_with_regions
from get_reads_positions import reads_positions
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
# matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
# matplotlib.rcParams["font.family"] = "sans-serif"
from scipy import stats
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})

sns.set_style("ticks")
import CTCF_TALL_modules_new
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"



def main():

    jurkat_df = pd.read_csv('f2_combined_HM_RPKM_csv/JURKAT_H3K27ac_combined.csv',index_col=0)
    cd4_df = pd.read_csv('f2_combined_HM_RPKM_csv/CD4_H3K27ac_combined.csv',index_col=0)

    union_feature_df = CTCF_TALL_modules_new.return_cancer_specific_combined_features('T-ALL')
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    gained_df = union_feature_df.loc[gained_df.index]
    const_df = CTCF_TALL_modules_new.return_constitutive_df()
    #print(gained_df);exit()
    
    ## Gained/Lost/Ctrl
    figname = 'T-ALL_H3K27ac_compr.pdf'
    fig = plt.figure(figsize = (2.8,2.8))
    i=0
    labels = ['Constitutive','T-ALL lost','T-ALL gained']
    colors = ['silver','blue','red']   
    for df_tmp in [const_df,lost_df,gained_df]:
        df1,df2 = jurkat_df.loc[df_tmp.index],cd4_df.loc[df_tmp.index]
        df_mean = 0.5*(df1.mean(axis=1)+df2.mean(axis=1))
        df_fc = np.log2(df1.mean(axis=1)+0.01)-np.log2(df2.mean(axis=1)+0.01)
        s,p = stats.ttest_ind(df1,df2,axis=1)
        a = plt.scatter(df_fc,-1*np.log10(p),label=labels[i],color=colors[i],s=20)
        #a = plt.scatter(df_mean,df_fc,label=labels[i],color=colors[i],s=6)
        i+=1
    plt.xlabel('H3K27ac ChIP-seq log$_2$(fold change)\nT-ALL over T-cell',fontsize=15)
    plt.ylabel('-log$_{{10}}$ $P$ value',fontsize=16)
    plt.xlim([-5.1,5.1])
    plt.ylim([-1,15])
#     plt.title('H3K27ac',fontsize=16)
    #plt.axes().tick_params(axis='y',direction='out', length=3, width=.8, colors='black')
    legend = plt.legend(frameon=False,fontsize=14,bbox_to_anchor=[.96,1],borderaxespad=0,labelspacing=.1,handletextpad=0.1,loc="upper left",markerscale=1.6) 
    legend.get_frame().set_facecolor("w")
    legend.get_frame().set_edgecolor("w")
#     sns.despine(offset=None, trim=False)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()
    

    ## Gained with dynamic H3K27ac vs. others
    figname = 'T-ALL_H3K27ac_compr_Gained_with_DynamicNotch.pdf'
    fig = plt.figure(figsize = (2.8,2.8))
    i=0
    labels = ['Constitutive','\nT-ALL gained with\ndynamic NOTCH1']
    colors = ['silver','red']  
    combined_df = pd.concat([pd.concat([const_df,lost_df]),gained_df]) 
#     gained_with_dynamic_notch = gained_df.loc[[553980, 465170, 304797, 529588, 91342, 22689, 119233, 403826, 89593, 349502, 539603, 219658, 266720, 358861, 370177, 320340, 562601, 275393, 222972]]
    gained_with_dynamic_notch = gained_df[gained_df['if_intra_domain_dynamic_notch']==1]
    for df_tmp in [combined_df,gained_with_dynamic_notch]:
        df1,df2 = jurkat_df.loc[df_tmp.index],cd4_df.loc[df_tmp.index]
        df_mean = 0.5*(df1.mean(axis=1)+df2.mean(axis=1))
        df_fc = np.log2(df1.mean(axis=1)+0.01)-np.log2(df2.mean(axis=1)+0.01)
        s,p = stats.ttest_ind(df1,df2,axis=1)
        a = plt.scatter(df_fc,-1*np.log10(p),label=labels[i],color=colors[i],s=20)
        #a = plt.scatter(df_mean,df_fc,label=labels[i],color=colors[i],s=6)
        i+=1
    plt.xlabel('H3K27ac ChIP-seq log$_2$(fold change)\nT-ALL over T-cell',fontsize=15)
    plt.ylabel('-log$_{{10}}$ $p$-value',fontsize=16)
    plt.xlim([-5.1,5.1])
    plt.ylim([-1,15])
#     plt.title('H3K27ac',fontsize=16)
    #plt.axes().tick_params(axis='y',direction='out', length=3, width=.8, colors='black')
    legend = plt.legend(frameon=False,fontsize=14,bbox_to_anchor=[.96,1],borderaxespad=0,labelspacing=-.5,handletextpad=0.1,loc="upper left",markerscale=1.6) 
    legend.get_frame().set_facecolor('w')
    legend.get_frame().set_edgecolor('w')
#     sns.despine(offset=None, trim=False)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
