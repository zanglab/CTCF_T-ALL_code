import os,sys,argparse,glob
import numpy as np
import pandas as pd
import find_overlap_keep_info_NOT_sep_strand_asimport
import re
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=18
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
import CTCF_TALL_modules_new

#>>> set(df['annotation'].values)
#{'Exon', 'Promoter', "5' UTR", 'Distal', "3' UTR", 'Intron'}


    
def sns_heatmap(df,title,figname):
    
    plt.figure(figsize=(2.6,2.6))
    g = sns.heatmap(data=df,cmap=plt.cm.GnBu,square=True,cbar_kws={"shrink":0.7},)
    cbar = g.collections[0].colorbar
    cbar.set_ticks([0,1])
    g.set_xticklabels([' '.join(i.split('_')) for i in df.columns],rotation=30,ha='right')
    g.set_yticklabels([' '.join(i.split('_')) for i in df.columns[::-1]],rotation=0)
    plt.title(title,fontsize=20)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,transparent=True)
    plt.close()



def main():
    
    outdir = 'f1_binding_overlap'
    os.makedirs(outdir,exist_ok=True)
    
    cancertypes=['T-ALL','AML','BRCA','CRC','LUAD','PRAD_TissueAdded']
    ticklabels=['T-ALL','AML','BRCA','CRC','LUAD','PRAD']
    

    gained_overlapped_matrix_raw = pd.DataFrame()
    lost_overlapped_matrix_raw = pd.DataFrame()
    gained_overlapped_matrix_jaccard = pd.DataFrame()
    lost_overlapped_matrix_jaccard = pd.DataFrame()   
    
    for cancertype1 in cancertypes:
        gained_df1,lost_df1 = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype1)
        for cancertype2 in cancertypes:
            gained_df2,lost_df2 = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype2)
            gained_overlapped = gained_df1.index.intersection(gained_df2.index)
            lost_overlapped = lost_df1.index.intersection(lost_df2.index)
            gained_common = gained_df1.index.union(gained_df2.index)
            lost_common = lost_df1.index.union(lost_df2.index)
            
            gained_overlapped_matrix_raw.loc[cancertype1,cancertype2] = len(gained_overlapped)
            lost_overlapped_matrix_raw.loc[cancertype1,cancertype2] = len(lost_overlapped)
            
            gained_overlapped_matrix_jaccard.loc[cancertype1,cancertype2] = len(gained_overlapped)/len(gained_common)
            lost_overlapped_matrix_jaccard.loc[cancertype1,cancertype2] = len(lost_overlapped)/len(lost_common)
    
    # raw 
    gained_overlapped_matrix_raw.to_csv('{}/panCancer_gained_CTCF_overlapped_raw.csv'.format(outdir))
    lost_overlapped_matrix_raw.to_csv('{}/panCancer_lost_CTCF_overlapped_raw.csv'.format(outdir))
    # jarccard
    gained_overlapped_matrix_jaccard.to_csv('{}/panCancer_gained_CTCF_overlapped_jaccard.csv'.format(outdir))
    lost_overlapped_matrix_jaccard.to_csv('{}/panCancer_lost_CTCF_overlapped_jaccard.csv'.format(outdir))
   

    # heat map of gained df
    gained_df_raw = pd.read_csv('{}/panCancer_gained_CTCF_overlapped_raw.csv'.format(outdir),index_col=0)
    gained_df = pd.read_csv('{}/panCancer_gained_CTCF_overlapped_jaccard.csv'.format(outdir),index_col=0)
    mask =  np.tri(gained_df.shape[1], k=-1)
    plt.figure(figsize=(3,3))
    data = gained_df_raw+np.sign(gained_df)*10
    data = data.clip(upper=110)
    g = sns.heatmap(data=data,mask=mask.T,cmap=plt.cm.Reds,square=True,cbar=False,cbar_kws={"shrink":0.62},vmax=200,vmin=1)
    g.set_yticklabels(ticklabels,rotation=0,ha='right',fontsize=16)
    g.set_xticklabels(ticklabels,rotation=30,ha='left',fontsize=16)
    plt.axes().tick_params(axis='x',direction='out', length=0, width=1, colors='black')
    plt.axes().tick_params(axis='y',direction='out', length=0, width=1, colors='black')
    g.xaxis.set_ticks_position('top') 
    for x in range(0,len(cancertypes)):
        for y in range(0,len(cancertypes)-x):
            plt.text(x + 0.05 , len(cancertypes)-y - 0.7, int(gained_df_raw.loc[gained_df.index[x],gained_df.columns[len(cancertypes)-1-y]]),fontsize=10,ha='left') #data[y,x] +0.05 , data[y,x] + 0.05
#     plt.title('Gained',ha='center',x=0.1,y=0.1,fontsize=20)
#     plt.axes().axhline(y=0, color='grey',linewidth=2)
#     plt.axes().axvline(x=0, color='grey',linewidth=2)
    plt.axes().invert_yaxis()
    plt.savefig('{}/panCancer_gained_CTCF_overlapped_jaccard.png'.format(outdir),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()
    
    
    # heat map of lost df
    lost_df_raw = pd.read_csv('{}/panCancer_lost_CTCF_overlapped_raw.csv'.format(outdir),index_col=0)
    lost_df = pd.read_csv('{}/panCancer_lost_CTCF_overlapped_jaccard.csv'.format(outdir),index_col=0)
    mask =  np.tri(lost_df.shape[1], k=-1)
    plt.figure(figsize=(2.8,2.8))
    data = lost_df_raw+np.sign(lost_df)*10
    data = data.clip(upper=40)
    g = sns.heatmap(data=data,mask=mask,cmap=plt.cm.Reds,yticklabels=True,square=True,cbar=False,cbar_kws={"shrink":0.7},vmax=80,vmin=1)
    g.set_xticklabels(ticklabels,rotation=30,ha='right',fontsize=16)
    g.set_yticklabels(ticklabels,rotation=0,ha='left',fontsize=16)
    plt.axes().tick_params(axis='x',direction='out', length=0, width=1, colors='black')
    plt.axes().tick_params(axis='y',direction='out', length=0, width=1, colors='black')
    g.yaxis.set_ticks_position('right') 
    for x in range(lost_df.shape[1]):
        for y in range(len(cancertypes)-1-x,lost_df.shape[1]):
            plt.text(x + 0.1 , len(cancertypes)-y - .7, int(lost_df_raw.loc[lost_df.index[x],lost_df.columns[len(cancertypes)-1-y]]),fontsize=10,ha='left') #data[y,x] +0.05 , data[y,x] + 0.05
    #cbar = g.collections[0].colorbar
    #cbar.set_ticks([0,4])
    #cbar.set_ticklabels([0,1])
#     plt.axes().axhline(y=4, color='grey',linewidth=2)
#     plt.axes().axvline(x=4, color='grey',linewidth=2)
#     plt.title('Lost',ha='center',x=0.9,y=4.7,fontsize=20)
    plt.axes().invert_yaxis()
    plt.savefig('{}/panCancer_lost_CTCF_overlapped_jaccard.png'.format(outdir),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()
    
    

           
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--infile1', action = 'store', type = str,dest = 'infile1', help = 'input file to be compared/separated', metavar = '<file>')
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
