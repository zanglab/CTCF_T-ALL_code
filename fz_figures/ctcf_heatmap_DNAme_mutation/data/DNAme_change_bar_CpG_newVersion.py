import os,sys,argparse
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=18
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.6)
sns.set_style("whitegrid", {'axes.grid' : False})
# from sklearn.metrics import r2_score
# from scipy import stats
# from scipy.stats import kde
# from scipy.interpolate import interpn
import re
import matplotlib.patches as mpatches
import scipy
import pandas as pd
from matplotlib import gridspec
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
import CTCF_TALL_modules_new





def plot_methylation_bar(df,dataType,outdir,flag):
    
    columns = ['differential_methy_EachSide150bp','DE_count_EachSide150bp']
    
    # ==== for those loci with #CpG < thre, treat as non-coveraged 
    diff_CpG_thre_cov = 3
    diff_CpG_thre_DMR = 3
    # ==== set diff. me. thre to show the bar
    avg_diff_me_thre = 20 
    
#     df = df[columns]
    df1 = df[df[columns[1]]>=diff_CpG_thre_cov]
#     df1[columns[0]][df[columns[1]]<diff_CpG_thre_DMR]=0
    df1 = df1.sort_values(by=columns[0],ascending=False)#;print(cancertype,flag,df)
    
    df2 = df[df[columns[1]] < diff_CpG_thre_cov]
    df = pd.concat([df1,df2])
    
    plt.figure(figsize = (2,5))
    norm = matplotlib.colors.Normalize(vmin=-100, vmax=100)
    color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=plt.cm.PiYG_r)
    
    for ii in np.arange(len(df.index)):
        # >=3 CpGs, with bar plot
        if df.loc[df.index[ii],columns[1]]>=diff_CpG_thre_DMR:
            plot_bar = df.loc[df.index[ii],columns[0]]#;print(plot_bar)
            if plot_bar>avg_diff_me_thre:
                g1 = plt.barh(ii,plot_bar,color='purple',height=1,lw=0)
            elif plot_bar<-1*avg_diff_me_thre:
                g1 = plt.barh(ii,plot_bar,color='green',height=1,lw=0)
        # with >=1 CpGs, but <3 CpGs, and plot with bg color
#         elif df.loc[df.index[ii],columns[1]]>=diff_CpG_thre_cov:
#             plot_bar = df.loc[df.index[ii],columns[0]]#;print(plot_bar)
#             if plot_bar>avg_diff_me_thre:
#                 g1 = plt.barh(ii,plot_bar,color='purple',height=1,lw=0)
#             elif plot_bar<-1*avg_diff_me_thre:
#                 g1 = plt.barh(ii,plot_bar,color='green',height=1,lw=0)
#             g2 = plt.barh(ii,100,color='lightgrey',height=1,lw=0,zorder=0)
#             g2 = plt.barh(ii,-100,color='lightgrey',height=1,lw=0,zorder=0)
        # un-covered regions   
        else:
            g2 = plt.barh(ii,-100,color='lightgrey',height=1,lw=0)
            g2 = plt.barh(ii,100,color='lightgrey',height=1,lw=0)
    
    
    plt.xlim([-100,100])
    plt.ylim([-.5,ii])
    plt.axes().text(100,df.shape[0]*1.3,'+',ha='center',fontsize=18)
    plt.axes().text(-100,df.shape[0]*1.3,'-',ha='center',fontsize=18)
    plt.axes().invert_yaxis()
    plt.axvline(0,c='k',lw=1.5)
    #plt.axvline(20,c='red')
    #plt.axvline(-20,c='red')
    plt.axes().set_yticks([])
    plt.ylabel('DNA methylation\nchanges',va='baseline',fontsize=18)
    #     plt.xlabel('+/-150bp')
    plt.axes().spines['left'].set_visible(False)
    plt.axes().spines['right'].set_visible(False)

    
    plt.savefig('{}/figs/{}_{}.pdf'.format(outdir,flag,dataType),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()
    return df
    


def main():

    outdir = 'plot_data_DNAme_changes_rank_by_methylation'
    os.makedirs(outdir,exist_ok=True)
    os.makedirs(outdir+os.sep+'figs',exist_ok=True)
    # here is to rank the cancertypes in the combined heatmap
#     cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']

    methylation_cancertype_matches = {\
'A549_cov5_vs_lung_cov5':'LUAD_2',\
'A549_cov5_vs_lung_cov5_WGBS':'LUAD_1',\
'CUTLL1_vs_A6010':'T-ALL_2',\
'HCT116_cov5_vs_sigmoid_colon_cov5':'CRC_1',\
'HCT116_GSM1465024_cov5_WGBS_vs_sigmoid_colon_cov5':'CRC_2',\
'HCT116_GSM257964x_cov5_WGBS_vs_sigmoid_colon_cov5':'CRC_3',\
'Jurkat_vs_A6010':'T-ALL_1',\
'LNCaP_cov5_vs_prostate_gland_cov5':'PRAD_TissueAdded',\
'MCF7_cov5_vs_breast_cov5':'BRCA',\
'PD31_cov5_vs_A6010':'T-ALL_patient1',\
'PD9_cov5_vs_A6010':'T-ALL_patient2'}


   
    for dataType in methylation_cancertype_matches.keys():
        cancertype=methylation_cancertype_matches[dataType]
        binding_cancertype = cancertype.split('_')[0]
        if binding_cancertype =="PRAD":
            binding_cancertype = 'PRAD_TissueAdded'
        gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding(binding_cancertype)
        methylation_info_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f8_DNA_methylation/f2_union_binding_methylation_compr/f3_union_binding_differential_DNA_methylation_append/{}.csv".format(dataType)
        print(os.path.isfile(methylation_info_file))
        methylation_info_df = pd.read_csv(methylation_info_file,index_col=0)
        df_plot = plot_methylation_bar(methylation_info_df.loc[gained_df.index],dataType,outdir,'gained')
        df_plot.to_csv(outdir+os.sep+'gained_{}.csv'.format(methylation_cancertype_matches[dataType]))
        
        df_plot = plot_methylation_bar(methylation_info_df.loc[lost_df.index],dataType,outdir,'lost')
        df_plot.to_csv(outdir+os.sep+'lost_{}.csv'.format(methylation_cancertype_matches[dataType]))
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

