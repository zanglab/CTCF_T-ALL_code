import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
matplotlib.use('Agg')
from matplotlib import gridspec
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=17
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
from matplotlib.colors import LinearSegmentedColormap
# import association_with_genes
# import association_with_regions
import re,bisect
# import CTCF_TALL_modules_new
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k',"axes.linewidth": .3})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]


def plot_mutation_changes(df,gs,loc):    
    ax = plt.subplot(gs[0,loc])
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
            g1 = ax.barh(ii,plot_bar,color='purple',height=1.5,edgecolor=None,linewidth=0)
        elif plot_bar<0:
            g1 = ax.barh(ii,plot_bar,color='green',height=1.5,edgecolor=None,linewidth=0)
    ax.set_xlim([xleft,xright])
#     plt.axes().set_xticks([xleft,xright])
#     plt.axes().set_xticklabels([xleft,xright])
    ax.set_ylim([-.5,ii])
#     plt.axes().text(xright,df.shape[0]*1.3,'+',ha='center',fontsize=20)
#     plt.axes().text(xleft,df.shape[0]*1.3,'-',ha='center',fontsize=20)
    ax.invert_yaxis()
    ax.axvline(0,c='k',lw=1.1)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel('$\Delta$(motif score)',va='baseline',fontsize=18)
#     ax.text(xright*1.2,df.shape[0],'{}'.format(df.shape[0]),fontsize=16,ha='left')
#     plt.xlabel('+/-150bp')
    ax.spines['left'].set_visible(True)
    ax.spines['right'].set_visible(True)




def sub_heatmap_plot(df,gs,loc,title,total,satuation,vmin,vmax,flag,cancertype):
    # plot each heatmap panel
    all_values = [i for col in df.columns for i in df[col].values]
    df = df.clip(upper=np.percentile(all_values,satuation))#;print(df)
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
        ax.set_title('CTCF\n{}'.format(title),fontsize=18) 
        if title=='T cell':
            ax.set_title('\n{}'.format(title),fontsize=18)    
        ax.set_ylabel('CTCF sites'.format(flag.capitalize()),fontsize=18)
        # only show the first/last xtick
        xp = g.get_xticks()
        ax.set_xticks([xp[0],xp[-1]])
        ax.set_xticklabels(['',''])
        ax.tick_params(axis='x',direction='out', length=0, width=1, colors='black')    
    
    # ==== show #bindings near the last heat map
    else:
        g=sns.heatmap(df,ax=ax,yticklabels=False,xticklabels=False,cbar=True,cmap=pal,vmin=cbarvmin,vmax = vmax,cbar_kws={"shrink": 0.5})
        ax.set_title('CTCF\n{}'.format(title),fontsize=18)   
        ax.set_ylabel('')
        cbar = g.collections[0].colorbar    
        cbar.set_ticks([cbarvmin,vmax])
        cbar.set_ticklabels([vmin,vmax])
        cbar.remove()
        



def plot_methylation_bar(df,gs,loc,cancertype,flag):

    ax = plt.subplot(gs[0,loc])

    columns = ['differential_methy_EachSide150bp','DE_count_EachSide150bp','UnionScored_count_EachSide150bp',]
    # ==== for those loci with #CpG < thre, treat as non-coveraged 
    diff_CpG_thre_cov = 3
    # ==== set diff. me. thre to show the bar
    avg_diff_me_thre = 20 
    
    # ==== confident DMR
    df1 = df[df[columns[1]]>=diff_CpG_thre_cov]
    df1 = df1.sort_values(by=columns[0],ascending=False)
    
    # ==== un-confident DMR (df2) and un-covered region (df3)
    df3 = df[df[columns[2]] ==0]
    df2 = df.loc[df.index.difference(df1.index.union(df3.index))]
    df = pd.concat([df1,df2,df3])
    bar_pos = df1.shape[0]+df2.shape[0]


    for ii in np.arange(len(df.index)):
        # >=3 CpGs, with bar plot
        if df.loc[df.index[ii],columns[1]]>=diff_CpG_thre_cov:
            plot_bar = df.loc[df.index[ii],columns[0]]
            if plot_bar>avg_diff_me_thre:
                g1 = ax.barh(ii,plot_bar,color='purple',height=1,lw=0)
            elif plot_bar<-1*avg_diff_me_thre:
                g1 = ax.barh(ii,plot_bar,color='green',height=1,lw=0)
        
        elif df.loc[df.index[ii],columns[2]] !=0 :
            g2 = ax.barh(ii,-100,color='lightgrey',height=1,lw=0)
            g2 = ax.barh(ii,100,color='lightgrey',height=1,lw=0)
        
        else:
            g2 = ax.barh(ii,-97,color='lightgrey',height=1,lw=0,zorder=9)
            g2 = ax.barh(ii,97,color='lightgrey',height=1,lw=0,zorder=9)
    

    ax.axhline(bar_pos-0.5,c='k',lw=1.1,zorder=10)
    ax.vlines(0,ymin=0-0.5,ymax=bar_pos-0.5,lw=1.1)

    ax.set_xlim([-100,100])
    ax.set_ylim([-.5,ii])
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_title('$\Delta$ DNA\n    methylation',fontsize=18,x=.75)
    
    # ==== add additional title
    add_title='{} {} sites (n={})'.format(cancertype.split('_')[0],flag,df.shape[0])
    ax.text(-700,-0.3*ii-.5,add_title,fontsize=18)    
    ax.hlines(y=-0.27*ii-.5,xmin=-680,xmax=180,clip_on=False,lw=1.1)

    ax.spines['left'].set_visible(True)
    ax.spines['right'].set_visible(True)
#     ax.text(115,df.shape[0]-1,'{:4}'.format(df.shape[0]),fontsize=15,ha='left')
    return df







###

# methylation_cancertype_matches = {\
# 'A549_cov5_vs_lung_cov5':'LUAD_2',\
# 'A549_cov5_vs_lung_cov5_WGBS':'LUAD_1',\
# 'CUTLL1_vs_A6010':'T-ALL_2',\
# 'HCT116_cov5_vs_sigmoid_colon_cov5':'CRC_1',\
# 'HCT116_GSM1465024_cov5_WGBS_vs_sigmoid_colon_cov5':'CRC_2',\
# 'HCT116_GSM257964x_cov5_WGBS_vs_sigmoid_colon_cov5':'CRC_3',\
# 'Jurkat_vs_A6010':'T-ALL_1',\
# 'LNCaP_cov5_vs_prostate_gland_cov5':'PRAD',\
# 'MCF7_cov5_vs_breast_cov5':'BRCA',\
# 'PD31_cov5_vs_A6010':'T-ALL_patient1',\
# 'PD9_cov5_vs_A6010':'T-ALL_patient2'}


cancertype_title_list = \
    {'T-ALL':{'methy_matches':"T-ALL_1",   'titles':['normal','cancer']},\
     'AML':	 {'methy_matches':"None",	 'titles':['normal','cancer']},\
     'BRCA': {'methy_matches':"BRCA",	 'titles':['normal','cancer']},\
     'CRC':	 {'methy_matches':"CRC_1",	 'titles':['normal','cancer']},\
     'LUAD': {'methy_matches':"LUAD_1",	 'titles':['normal','cancer']},\
     'PRAD_TissueAdded': {'methy_matches':"PRAD_TissueAdded",	 'titles':['normal','cancer']}}


outdir = 'pancancer_fig_DNAme'
os.makedirs(outdir,exist_ok=True)
os.makedirs(outdir+os.sep+'csv',exist_ok=True)

cancertypes = ["T-ALL","CRC","BRCA","LUAD","PRAD_TissueAdded"]
bindingtypes = ["gained","lost"]
satuation,vmin,vmax = 98,0,100

# cancertypes = ["BRCA"]
# bindingtypes = ["gained","lost"]

for bindingtype in bindingtypes:
    for cancertype in cancertypes:
  
        normal_file = "data/plot_data_CTCF_heatmap_rerank_by_methylation/{}_{}_normal.csv".format(bindingtype,cancertype)
        cancer_file = "data/plot_data_CTCF_heatmap_rerank_by_methylation/{}_{}_cancer.csv".format(bindingtype,cancertype)
        mutation_file = "data/plot_data_mutation/{}_{}.csv".format(bindingtype,cancertype)     
        methylation_file = "data/plot_data_DNAme_changes_rank_by_methylation/{}_{}.csv".format(bindingtype,cancertype_title_list[cancertype]['methy_matches'])     

        
        #### plot
        fig = plt.figure(figsize = (3,3))
        total = 4
        width_ratio = [1.1,1.1,.0,.8]
        gs = gridspec.GridSpec(1,total,width_ratios=width_ratio,wspace=0) 
        

        #### subplot of methylation
        loc = 3
        df = pd.read_csv(methylation_file,index_col=0)
        re_index_df = plot_methylation_bar(df,gs,loc,cancertype,bindingtype)
        re_index_df.to_csv(outdir+'/csv/{}_{}_methylation.csv'.format(bindingtype,cancertype))
        

        #### subplot of normal CTCF ChIP 
        loc = 0
        df = pd.read_csv(normal_file,index_col=0)
        df = df.reindex(re_index_df.index)
        df.to_csv(outdir+'/csv/{}_{}_CTCF_normal.csv'.format(bindingtype,cancertype))
        sub_heatmap_plot(df,gs,loc,cancertype_title_list[cancertype]['titles'][loc],total,satuation,vmin,vmax,bindingtype,cancertype)


        #### subplot of cancer CTCF ChIP 
        loc = 1
        df = pd.read_csv(cancer_file,index_col=0)
        df = df.reindex(re_index_df.index)
        df.to_csv(outdir+'/csv/{}_{}_CTCF_cancer.csv'.format(bindingtype,cancertype))
        sub_heatmap_plot(df,gs,loc,cancertype_title_list[cancertype]['titles'][loc],total,satuation,vmin,vmax,bindingtype,cancertype)

        
#         loc = 5
#         df = pd.read_csv(mutation_file,index_col=0)        
#         plot_mutation_changes(df,gs,loc)

        plt.savefig(outdir+'/{}_{}.pdf'.format(bindingtype,cancertype),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
        plt.close()
        
#         exit()
        

        




    
