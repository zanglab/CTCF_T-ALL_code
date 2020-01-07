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
sns.set(font_scale=2)
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


        
def rank_index_by_methylation_figs(df,cancertype,flag):
    
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype)
    columns = ['differential_methy_EachSide150bp','DE_count_EachSide150bp']
    if columns[0] in df.columns:
        if flag=='gained':
            rank_index = gained_df.sort_values(by=columns[0],ascending=False).index
        elif flag=='lost':
            rank_index = lost_df.sort_values(by=columns[0],ascending=False).index
    else:
        if flag=='gained':
            rank_index = gained_df.sort_values(by='cancer_vs_other_stats',ascending=False).index
        elif flag=='lost':
            rank_index = lost_df.sort_values(by='cancer_vs_other_stats',ascending=True).index

    df = df.loc[rank_index]
    return df
        

def plot_mutation_changes(df,cancertype,outdir,flag):    
    plt.figure(figsize = (1.3,2.8))
    #norm = matplotlib.colors.Normalize(vmin=-100, vmax=100)
    #color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=plt.cm.PiYG_r)
    df['mutation_change']=0
    df.loc[df[df['strand']=='+'].index,'mutation_change']=df['alt_seq_score']-df['seq_score']
    df.loc[df[df['strand']=='-'].index,'mutation_change']=df['rev_alt_seq_score']-df['rev_seq_score']
    df['mutation_change'] = df['mutation_change'].fillna(0)
    # to give the column name for figs
    columns=[['mutation_change']]
    
    xleft,xright = -9,9
    for ii in np.arange(len(df.index)):
        plot_bar = df.loc[df.index[ii],columns[0]].values[0]#;print(plot_bar)
        if plot_bar>0:
            g1 = plt.barh(ii,plot_bar,color='purple',height=1.5,edgecolor=None,linewidth=0)
        elif plot_bar<0:
            g1 = plt.barh(ii,plot_bar,color='green',height=1.5,edgecolor=None,linewidth=0)
    plt.xlim([xleft,xright])
    plt.axes().set_xticks([xleft,xright])
    plt.axes().set_xticklabels([xleft,xright])
    plt.ylim([-.5,ii])
    plt.axes().text(xright,df.shape[0]*1.3,'+',ha='center',fontsize=20)
    plt.axes().text(xleft,df.shape[0]*1.3,'-',ha='center',fontsize=20)
    plt.axes().invert_yaxis()
    plt.axvline(0,c='k',lw=1.5)
#     plt.axvline(20,c='r')
#     plt.axvline(-20,c='r')
    plt.axes().set_yticks([])
    plt.ylabel('$\Delta$(motif score)',va='baseline',fontsize=20)
#     plt.xlabel('+/-150bp')
    plt.axes().spines['left'].set_visible(False)
    plt.axes().spines['right'].set_visible(False)
    plt.savefig('{}/{}_{}.pdf'.format(outdir,flag,cancertype),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()
    df.to_csv('{}/{}_{}.csv'.format(outdir,flag,cancertype))
    return df
    

def main():

    outdir = 'f7_mutation_change_bar'
    os.makedirs(outdir,exist_ok=True)
    # here is to rank the cancertypes in the combined heatmap
    cancertype_mutation_matchness = {'BRCA':'BRCA','CRC':'COAD','LUAD':'LUAD','PRAD':'PRAD','AML':'AML','PRAD_TissueAdded':'PRAD'}
    cancertypes=['BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']
    for cancertype in cancertypes:
        for binding_type in ['gained','lost']:
            mutation_change_file = 'f6_mutation_matrix_score/{}_{}.csv'.format(cancertype,binding_type)
            mutation_df = pd.read_csv(mutation_change_file,index_col=0)
            mutation_df = rank_index_by_methylation_figs(mutation_df,cancertype,binding_type)
            mutation_df = plot_mutation_changes(mutation_df,cancertype,outdir,binding_type)
                






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

