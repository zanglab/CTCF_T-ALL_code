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
from sklearn.metrics import r2_score
from scipy import stats
from scipy.stats import kde
from scipy.interpolate import interpn
import re
import matplotlib.patches as mpatches
import scipy
import pandas as pd
from matplotlib import gridspec
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
import CTCF_TALL_modules_new





def plot_methylation_bar(df,dataType,outdir,flag):
    
    columns = ['differential_methy_EachSide150bp','DE_count_EachSide150bp']
    df = df[columns]
    df[columns[0]][df[columns[1]]<5]=0
    df = df.sort_values(by=columns[0],ascending=False)#;print(cancertype,flag,df)
    
    plt.figure(figsize = (1.3,2.8))
    norm = matplotlib.colors.Normalize(vmin=-100, vmax=100)
    color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=plt.cm.PiYG_r)
    for ii in np.arange(len(df.index)):
        plot_bar = df.loc[df.index[ii],columns[0]]#;print(plot_bar)
    #    plot_bar_cov = df.loc[df.index[ii],'DE_count_EachSide150bp_cov5']
        if plot_bar>0:
            g1 = plt.barh(ii,plot_bar,color='purple',height=1,lw=0)
            #g2 = plt.barh(ii,100,color='antiquewhite',height=1,zorder=0)
        elif plot_bar<0:
            g1 = plt.barh(ii,plot_bar,color='green',height=1,lw=0)
            #g2 = plt.barh(ii,-100,color='antiquewhite',height=1,zorder=0)
    plt.xlim([-100,100])
    plt.ylim([-.5,ii])
    plt.axes().text(100,df.shape[0]*1.3,'+',ha='center',fontsize=18)
    plt.axes().text(-100,df.shape[0]*1.3,'-',ha='center',fontsize=18)
    plt.axes().invert_yaxis()
    plt.axvline(0,c='k',lw=1.5)
    plt.axvline(20,c='red')
    plt.axvline(-20,c='red')
    plt.axes().set_yticks([])
    plt.ylabel('DNA methylation\nchanges',va='baseline',fontsize=18)
#     plt.xlabel('+/-150bp')
    plt.axes().spines['left'].set_visible(False)
    plt.axes().spines['right'].set_visible(False)
    
    plt.savefig('{}/{}_{}.pdf'.format(outdir,flag,dataType),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()
    


def main():

    outdir = 'f6_DNAme_change_bar_CpG_version_sep_cancertype'
    os.makedirs(outdir,exist_ok=True)
    # here is to rank the cancertypes in the combined heatmap
    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']

    methylation_cancertype_matches = {'A549_cov5_vs_lung_cov5':'LUAD',\
'A549_cov5_vs_lung_cov5_WGBS':'LUAD',\
'CUTLL1_vs_A6010':'T-ALL',\
'HCT116_cov5_vs_sigmoid_colon_cov5':'CRC',\
'HCT116_GSM1465024_cov5_WGBS_vs_sigmoid_colon_cov5':'CRC',\
'HCT116_GSM257964x_cov5_WGBS_vs_sigmoid_colon_cov5':'CRC',\
'Jurkat_vs_A6010':'T-ALL',\
'MCF7_cov5_vs_breast_cov5':'BRCA',\
'PD31_cov5_vs_A6010':'T-ALL',\
'PD9_cov5_vs_A6010':'T-ALL'}
    
    for dataType in methylation_cancertype_matches.keys():
        gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding(methylation_cancertype_matches[dataType])
        methylation_info_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f8_DNA_methylation/f2_union_binding_methylation_compr/f2_union_binding_differential_DNA_methylation/{}.csv".format(dataType)
        print(os.path.isfile(methylation_info_file))
        methylation_info_df = pd.read_csv(methylation_info_file,index_col=0)
        plot_methylation_bar(methylation_info_df.loc[gained_df.index],dataType,outdir,'gained')
        plot_methylation_bar(methylation_info_df.loc[lost_df.index],dataType,outdir,'lost')
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

