import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import association_with_genes
import association_with_regions
import re,bisect
import CTCF_TALL_modules_new

import matplotlib
matplotlib.use('Agg')
from matplotlib import gridspec
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set_style("whitegrid", {'axes.grid' : False,'grid.color': 'grey'})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'w'})
from matplotlib.colors import LinearSegmentedColormap
#plus = re.compile('\+')
#minus = re.compile('\-')
matplotlib.rcParams["font.sans-serif"] = ["Arial"]

cancer_type_matched = {'T-ALL':['CD4','JURKAT'],'BRCA':['Breast_normal','BRCA'],'LUAD':['Lung_normal','LUAD',],'CRC':['Colon_normal','CRC'],'AML':['Erythroid','AML']}
ctcf_cancer_type_matched = {'T-ALL':['T-ALL_normal','T-ALL_cancer'],\
    'AML':['AML_normal','AML_cancer'],\
    'BRCA':['BRCA_normal','BRCA_cancer'],\
    'CRC':['CRC_normal','CRC_cancer'],\
    'LUAD':['LUAD_normal','LUAD_cancer']}
    
cancer_type_matched_title = {'T-ALL':['T-cell','T-ALL'],'BRCA':['Breast \nTissue','BRCA',],'LUAD':['Lung \nTissue','LUAD',],'CRC':['Colon \nTissue','CRC',],'AML':['Erythroid','AML']}
cancer_type_matched_title = {'T-ALL':['normal','cancer'],'BRCA':['normal','cancer',],'LUAD':['normal','cancer',],'CRC':['normal','cancer',],'AML':['normal','cancer']}
    
def window_cumulative(df,half_window=7,step=1):
    smooth_df_columns = np.arange(0,len(df.columns),step)
    smooth_df = pd.DataFrame(index=df.index,columns=smooth_df_columns)#;print(smooth_df.columns)
    for col in smooth_df_columns:
        window_left = max(col-half_window,0)
        window_right = min(col+half_window,len(df.columns)-1)
        smooth_df.loc[:,col] = df.iloc[:,window_left:window_right+1].sum(axis=1)    
    #print(df,smooth_df)
    return smooth_df   

def signal_centered(df):
    center_position = int(df.shape[1]/2)
    for row in df.index:
        vals = df.loc[row]
        max_index = list(vals).index(vals.max())
        # move max to center
        if max_index<center_position:
            df.loc[row] = np.append(np.zeros(center_position-max_index),vals)[:df.shape[1]]
        elif max_index>center_position:
            df.loc[row] = np.append(vals[max_index-center_position:],np.zeros(max_index-center_position))
    return df

def window_smooth(df,half_window):
    smooth_df = pd.DataFrame(index=df.index,columns=df.columns)
    for col in np.arange(len(df.columns)):
        window_left = max(col-half_window,0)
        window_right = min(col+half_window,len(df.columns)-1)
        smooth_df.iloc[:,col] = df.iloc[:,window_left:window_right+1].mean(axis=1)
    return smooth_df   
    print(df,smooth_df);exit()
    return df


def return_CTCF_state_df(cancertype,cancer_subname,ctcf_subname,flag,hm):
    # return the df for heatmap
    binding_dir = 'f1_gained_lost_HM_signals'
    binding_file = binding_dir+os.sep+'{}_{}_{}.csv'.format(flag,cancer_subname,hm)
    if hm=='CTCF':
        binding_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f8_panCancer_CTCF_binding_pattern/f3_cancerType_CTCF_heatmap/f1_gained_lost_CTCF_binding'
        binding_file = binding_dir+os.sep+'{}_{}.csv'.format(flag,ctcf_subname)
    
    with open(binding_file) as binding_inf:
        binding_df = pd.read_csv(binding_inf,index_col=0)
    
#     ## rank by K27 sig
    loc_index_flag=0 if flag=='lost' else 1
    rank_hm = 'H3K27ac'
    index_binding_dir = 'f1_gained_lost_HM_signals'
    index_binding_file = index_binding_dir+os.sep+'{}_{}_{}.csv'.format(flag,cancer_type_matched[cancertype][loc_index_flag],rank_hm);print(index_binding_file)
    with open(index_binding_file) as index_binding_inf:
        index_binding_df = pd.read_csv(index_binding_inf,index_col=0)
    rank_index = index_binding_df.mean(axis=1).sort_values(ascending=False).index

#     gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype)
#     columns = ['differential_methy_EachSide150bp','DE_count_EachSide150bp']
#     if columns[0] in gained_df.columns:
#         if flag=='gained':
#             rank_index = gained_df.sort_values(by=columns[0],ascending=False).index
#         elif flag=='lost':
#             rank_index = lost_df.sort_values(by=columns[0],ascending=False).index
#     else:
#         if flag=='gained':
#             rank_index = gained_df.sort_values(by='cancer_vs_other_stats',ascending=False).index
#         elif flag=='lost':
#             rank_index = lost_df.sort_values(by='cancer_vs_other_stats',ascending=True).index
 
    binding_df = binding_df.loc[rank_index]
    binding_df = window_cumulative(binding_df)
    if hm=='CTCF':
        binding_df = signal_centered(binding_df)
    return binding_df


def sub_heatmap_plot(df,gs,title,loc,total,satuation,vmin,vmax,flag,cancertype,hm):
    # plot each heatmap panel
    all_values = [i for col in df.columns for i in df[col].values]
    df = df.clip(upper=np.percentile(all_values,satuation))
    ax = plt.subplot(gs[0,loc])
    #pal = sns.light_palette(cancertype_heatmap_color(cancertype),as_cmap=True)
    pal = sns.light_palette('blue',as_cmap=True)
    if hm=='CTCF':
        pal = sns.light_palette('red',as_cmap=True)
    cbarvmin=0
    if loc%3==2:
        ax.set_axis_off()
    elif loc==0:   
        g=sns.heatmap(df,ax=ax,yticklabels=False,xticklabels=True,cbar=True,cmap=pal,vmin=cbarvmin,vmax = vmax,cbar_kws={"shrink": 0.5})
        ax.set_ylabel('{} ChIP-seq'.format(hm),fontsize=13)
        xp = g.get_xticks()#;print(xp)
        ax.set_xticks([xp[0],xp[-1]])
        ax.set_xticklabels(['-1kb','1kb'],rotation=30,fontsize=13)
        ax.set_title('{}'.format(title),fontsize=14)    
        ax.tick_params(axis='x',direction='out', length=0, width=1, colors='black')    
        cbar = g.collections[0].colorbar
        cbar.set_clim(vmax*.15,vmax)
        cbar.remove()
    elif loc==total-1:
        g=sns.heatmap(df,ax=ax,yticklabels=False,xticklabels=False,cbar=True,cmap=pal,vmin=cbarvmin,vmax = vmax,cbar_kws={"shrink": 0.5})
        ax.set_title('{}'.format(title),fontsize=15)
        ax.set_ylabel('')
        ax.tick_params(axis='y',direction='out', length=0, width=1, colors='black')    
        cbar = g.collections[0].colorbar
        cbar.set_clim(vmax*.15,vmax)
        cbar.set_ticks([cbarvmin,vmax])
        cbar.set_ticklabels([vmin,vmax]) 
        cbar.ax.set_position([.9,0.35,.8,.5]) 
        cbar.ax.tick_params(axis='y',direction='out', length=0, width=1, colors='black')            
           
    else:
        g=sns.heatmap(df,ax=ax,yticklabels=False,xticklabels=False,cbar=True,cmap=pal,vmin=cbarvmin,vmax = vmax)
        ax.set_title('{}'.format(title),fontsize=15)
        ax.set_ylabel('')
        ax.tick_params(axis='y',direction='out', length=0, width=1, colors='black')    
        cbar = g.collections[0].colorbar
        cbar.set_clim(vmax*.15,vmax)
        cbar.remove()
        
   
    if hm == 'H3K27me3' and loc%3==1:
        ax.text(210,df.shape[0]*1.15,'{}'.format(df.shape[0]),fontsize=12,ha='left')
   

    if loc%3==1:
        # ==== add additional title
        ax.text(-100,-0.22*df.shape[0]-.3,cancertype,fontsize=15)    
        ax.hlines(y=-0.2*df.shape[0]-.3,xmin=-230,xmax=200,clip_on=False,lw=1.1)


def combine_plot(cancertypes,hm,flag,outdir,satuation,vmin,vmax):
    # for each cancer type, return cancer/normal df for heatmap
    total_len=sum([len(cancer_type_matched[cancertype]) for cancertype in cancertypes]) +len(cancertypes)-1
    #print(total_len);exit()
    fig = plt.figure(figsize = (0.8*total_len,1.8))
    # to make sure the last figure have more space for color bar
    width_ratio = np.append([1]*(total_len-1),[1.2])
    width_ratio = [1,1,.35,1,1,.35,1,1,.35,1,1,]
    gs = gridspec.GridSpec(1,total_len,width_ratios=width_ratio,wspace=-0.1) 
    i=0
    for cancertype in cancertypes:
        for loc_index in np.arange(len(cancer_type_matched[cancertype])):
            if i%3==2:
                sub_heatmap_plot(loc_df,gs,cancer_title,i,total_len,satuation,vmin,vmax,flag,cancertype,hm)
                i+=1
            cancer_subname = cancer_type_matched[cancertype][loc_index]
            ctcf_subname = ctcf_cancer_type_matched[cancertype][loc_index]
            cancer_title = cancer_type_matched_title[cancertype][loc_index]
            loc_df = return_CTCF_state_df(cancertype,cancer_subname,ctcf_subname,flag,hm)
            sub_heatmap_plot(loc_df,gs,cancer_title,i,total_len,satuation,vmin,vmax,flag,cancertype,hm)
            i+=1 
    plt.savefig(outdir+os.sep+'{}_{}_combi.png'.format(flag,hm),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close();


def main():

    outdir = 'f3_HM_heatmap_combined_PanCancer'
    os.makedirs(outdir,exist_ok=True)
    
    cancertypes = ['T-ALL','BRCA','CRC','LUAD']
    vmin=0
    hms = ['H3K27ac','H3K27me3','H3K4me1','CTCF']
    vmax_hm = [200,100,100,120]
    
    for ii in [0,1,2,3]:
        for satuation in [98]:
            for vmax in [vmax_hm[ii]]:
                combine_plot(cancertypes,hms[ii],'gained',outdir,satuation,vmin,vmax)
                combine_plot(cancertypes,hms[ii],'lost',outdir,satuation,vmin,vmax)
#                 exit()
            





    
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

