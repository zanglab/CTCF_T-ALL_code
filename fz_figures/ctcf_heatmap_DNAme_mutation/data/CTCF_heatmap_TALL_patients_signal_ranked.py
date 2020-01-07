import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
matplotlib.use('Agg')
from matplotlib import gridspec
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.6)
sns.set_style("whitegrid", {'axes.grid' : False})
from matplotlib.colors import LinearSegmentedColormap
import association_with_genes
import association_with_regions
import re,bisect
import CTCF_TALL_modules_new
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'w'})


cancertype_prenema_title_list = \
    {'T-ALL_p1':{'prenames':['T-ALL_normal','PD9'],'titles':['T cell','T-ALL\npatient1']},\
     'T-ALL_p2':{'prenames':['T-ALL_normal','PD31'],'titles':['T cell','T-ALL\npatient2']}}    

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



def return_CTCF_state_df(cancertype,prename,flag):
    # return the df for heatmap
    binding_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f8_panCancer_CTCF_binding_pattern/f3_cancerType_CTCF_heatmap/f1_gained_lost_CTCF_binding'
    binding_file = binding_dir+os.sep+'{}_{}.csv'.format(flag,prename)
    with open(binding_file) as binding_inf:
        binding_df = pd.read_csv(binding_inf,index_col=0)
    # ==== rank by patient signal
    patient_binding = pd.read_csv(binding_dir+os.sep+'{}_{}.csv'.format(flag,cancertype_prenema_title_list[cancertype]['prenames'][1]),index_col=0)
    if flag=='gained':
        rank_index = patient_binding.mean(axis=1).sort_values(ascending=False).index
    if flag=='lost':
        rank_index = patient_binding.mean(axis=1).sort_values(ascending=True).index
    binding_df = binding_df.loc[rank_index]
    binding_df = window_cumulative(binding_df)
    binding_df = signal_centered(binding_df)
    return binding_df



def sub_heatmap_plot(df,gs,loc,title,total,satuation,vmin,vmax,flag,cancertype):
    # plot each heatmap panel
    all_values = [i for col in df.columns for i in df[col].values]
    df = df.clip(upper=np.percentile(all_values,satuation))
    ax = plt.subplot(gs[0,loc])
    pal = sns.light_palette('red',as_cmap=True)
    cbarvmin=12
    # ==== show the y-label and x-ticklabel of the first heat map
    if loc==0:   
        g=sns.heatmap(df,ax=ax,yticklabels=False,xticklabels=True,cbar=True,cmap=pal,vmin=cbarvmin,vmax = vmax,cbar_kws={"shrink": 0.5})
        cbar = g.collections[0].colorbar    
        cbar.set_ticks([cbarvmin,vmax])
        cbar.set_ticklabels([vmin,vmax])
        cbar.remove()
        ax.set_title('{}'.format(title),fontsize=18)   
        ax.set_ylabel('CTCF ChIP-seq'.format(flag.capitalize()),fontsize=18)
        # only show the first/last xtick
        xp = g.get_xticks()
        ax.set_xticks([xp[0],xp[-1]])
        ax.set_xticklabels(['-1kb','1kb'],rotation=30,fontsize=16)
        ax.tick_params(axis='x',direction='out', length=0, width=1, colors='black')    
    
    # ==== show #bindings near the last heat map
    else:
        g=sns.heatmap(df,ax=ax,yticklabels=False,xticklabels=False,cbar=True,cmap=pal,vmin=cbarvmin,vmax = vmax,cbar_kws={"shrink": 0.5})
        cbar = g.collections[0].colorbar    
        cbar.set_ticks([cbarvmin,vmax])
        cbar.set_ticklabels([vmin,vmax])
        cbar.remove()
        ax.set_title('{}'.format(title),fontsize=18)   
        ax.set_ylabel('')
        if loc==total-1:
            ax.text(220,df.shape[0],'{}'.format(df.shape[0]),fontsize=17,ha='left')
        


def combine_plot(cancertype,flag,outdir,satuation,vmin,vmax):
    # for each cancer type, return cancer/normal df for heatmap
    prenames = cancertype_prenema_title_list[cancertype]['prenames']
    titles = cancertype_prenema_title_list[cancertype]['titles']
    total_len = len(prenames)
    fig = plt.figure(figsize = (3.2,2.7))
    width_ratio = [1,1]
    gs = gridspec.GridSpec(1,total_len,width_ratios=width_ratio,wspace=0.0) 
    for pos_id in np.arange(total_len):
        loc_df = return_CTCF_state_df(cancertype,prenames[pos_id],flag)
        loc_df.to_csv(outdir+os.sep+'{}_{}_{}.csv'.format(flag,cancertype,prenames[pos_id]))

        sub_heatmap_plot(loc_df,gs,pos_id,titles[pos_id],total_len,satuation,vmin,vmax,flag,cancertype)
    plt.savefig(outdir+os.sep+'/figs/{}_CTCF_{}.png'.format(flag,cancertype),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close();



def main():

    outdir = 'plot_data_CTCF_heatmap_rerank_by_methylation_patients'
    os.makedirs(outdir,exist_ok=True)
    os.makedirs(outdir+os.sep+'figs',exist_ok=True)

    cancertypes=['T-ALL_p1','T-ALL_p2']
    for cancertype in cancertypes:
        vmin=0
        for satuation in [98]:
            for vmax in [120]:
                # pass
                combine_plot(cancertype,'gained',outdir,satuation,vmin,vmax)#;exit()
                combine_plot(cancertype,'lost',outdir,satuation,vmin,vmax)


    # ==== heat map of color bar
#     vmax = 120
#     fig,ax = plt.subplots(figsize=(1.,.2))
#     norm = matplotlib.colors.Normalize(vmin=0, vmax=vmax)
#     pal = sns.light_palette('red',as_cmap=True)
#     cb = matplotlib.colorbar.ColorbarBase(ax,cmap=pal,norm=norm,orientation='horizontal')
# #     cb.set_label('change level',rotation=90,labelpad=-110)
#     cb.set_ticks([0,vmax])
#     cb.set_clim(vmax*.15,vmax)
#     ax.tick_params(axis='x',direction='out', length=0, width=1, colors='black')    
#     plt.savefig(outdir+os.sep+'colorbar.png',bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
#     plt.close()
# 
#     fig,ax = plt.subplots(figsize=(.25,1))
#     norm = matplotlib.colors.Normalize(vmin=0, vmax=vmax)
#     pal = sns.light_palette('red',as_cmap=True)
#     cb = matplotlib.colorbar.ColorbarBase(ax,cmap=pal,norm=norm,orientation='vertical')
# #     cb.set_label('change level',rotation=90,labelpad=-110)
#     cb.set_ticks([0,vmax])
#     cb.set_clim(vmax*.15,vmax)
#     plt.savefig(outdir+os.sep+'colorbar2.png',bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
#     plt.close()
# 




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

