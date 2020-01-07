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
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
from matplotlib.colors import LinearSegmentedColormap
# import association_with_genes
# import association_with_regions
import re,bisect
# import CTCF_TALL_modules_new
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k',"axes.linewidth": .3})


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
        ax.set_title('{}'.format(title),fontsize=18) 
        if title=='T cell':
            ax.set_title('\n{}'.format(title),fontsize=18)    
        ax.set_ylabel('CTCF ChIP-seq'.format(flag.capitalize()),fontsize=18,va='baseline')
        # only show the first/last xtick
        xp = g.get_xticks()
        ax.set_xticks([xp[0],xp[-1]])
        if cancertype.startswith('PRAD'):
            ax.set_xticklabels(['-1kb','1kb'],rotation=30,fontsize=16)
        else:
            ax.set_xticklabels(['',''])
        ax.tick_params(axis='x',direction='out', length=0, width=1, colors='black')    
    
    # ==== show #bindings near the last heat map
    else:
        g=sns.heatmap(df,ax=ax,yticklabels=False,xticklabels=False,cbar=True,cmap=pal,vmin=cbarvmin,vmax = vmax,cbar_kws={"shrink": 0.5})
        ax.set_title('{}'.format(title),fontsize=18)   
        ax.set_ylabel('')
        cbar = g.collections[0].colorbar    
        cbar.set_ticks([cbarvmin,vmax])
        cbar.set_ticklabels([vmin,vmax])
        cbar.remove()
        if loc==total-1:
            ax.text(220,df.shape[0],'{}'.format(df.shape[0]),fontsize=17,ha='left')
        



def plot_methylation_bar(df,gs,loc):
    
    columns = ['differential_methy_EachSide150bp','DE_count_EachSide150bp']
    df = df[columns]
    df[columns[0]][df[columns[1]]<5]=0
    df = df.sort_values(by=columns[0],ascending=False)#;print(cancertype,flag,df)
    
    ax = plt.subplot(gs[0,loc])
    norm = matplotlib.colors.Normalize(vmin=-100, vmax=100)
    color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=plt.cm.PiYG_r)
    for ii in np.arange(len(df.index)):
        plot_bar = df.loc[df.index[ii],columns[0]]#;print(plot_bar)
        if plot_bar>0:
            g1 = ax.barh(ii,plot_bar,color='purple',height=1,lw=0)
        elif plot_bar<0:
            g1 = ax.barh(ii,plot_bar,color='green',height=1,lw=0)
    ax.set_xlim([-100,100])
    ax.set_ylim([-.5,ii])
#     plt.axes().text(100,df.shape[0]*1.3,'+',ha='center',fontsize=18)
#     plt.axes().text(-100,df.shape[0]*1.3,'-',ha='center',fontsize=18)
    ax.invert_yaxis()
    ax.axvline(0,c='k',lw=1.1)
#     plt.axvline(20,c='red')
#     plt.axvline(-20,c='red')
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel('$\Delta$(DNA methylation)',va='baseline',fontsize=18)
    ax.spines['left'].set_visible(True)
    ax.spines['right'].set_visible(True)




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





###

cancertype_title_list = \
    {'T-ALL':{'prenames':['T-ALL_normal','T-ALL_cancer'],'titles':['T cell','T-ALL']},\
     'AML':	 {'prenames':['AML_normal','AML_cancer'],	 'titles':['Erythroid ','AML']},\
     'BRCA': {'prenames':['BRCA_normal','BRCA_cancer'],	 'titles':['BRCA\ntissue','BRCA']},\
     'CRC':	 {'prenames':['CRC_normal','CRC_cancer'],	 'titles':['CRC\ntissue','CRC']},\
     'LUAD': {'prenames':['LUAD_normal','LUAD_cancer'],	 'titles':['LUAD\ntissue','LUAD']},\
     'PRAD': {'prenames':['PRAD_normal','PRAD_cancer'],	 'titles':['PRAD\ntissue','PRAD']},\
     'PRAD_TissueAdded':{'prenames':['PRAD_TissueAdded_normal','PRAD_TissueAdded_cancer'],'titles':['PRAD\ntissue','PRAD']}}


outdir = 'pancancer_fig'
os.makedirs(outdir,exist_ok=True)

cancertypes = ["T-ALL","CRC","BRCA","LUAD"]
bindingtypes = ["gained","lost"]
satuation,vmin,vmax = 98,0,100


for bindingtype in bindingtypes:
    for cancertype in cancertypes:
  
        normal_file = "data/plot_data_CTCF_heatmap_rerank_by_methylation/{}_{}_normal.csv".format(bindingtype,cancertype)
        cancer_file = "data/plot_data_CTCF_heatmap_rerank_by_methylation/{}_{}_cancer.csv".format(bindingtype,cancertype)
        methylation_file = "data/plot_data_DNAme_changes_rank_by_methylation/{}_{}.csv".format(bindingtype,cancertype)     
        mutation_file = "data/plot_data_mutation/{}_{}.csv".format(bindingtype,cancertype)     
        
        #### plot
        fig = plt.figure(figsize = (3.9,3))
        total = 6
        width_ratio = [1,1,.35,.7,.5,.7]
        gs = gridspec.GridSpec(1,total,width_ratios=width_ratio,wspace=0) 
        
        #### subplot of normal CTCF ChIP 
        loc = 0
        df = pd.read_csv(normal_file,index_col=0)
        sub_heatmap_plot(df,gs,loc,cancertype_title_list[cancertype]['titles'][loc],total,satuation,vmin,vmax,bindingtype,cancertype)


        #### subplot of cancer CTCF ChIP 
        loc = 1
        df = pd.read_csv(cancer_file,index_col=0)
        sub_heatmap_plot(df,gs,loc,cancertype_title_list[cancertype]['titles'][loc],total,satuation,vmin,vmax,bindingtype,cancertype)

        #### subplot of methylation
        loc = 3
        df = pd.read_csv(methylation_file,index_col=0)
        plot_methylation_bar(df,gs,loc)
        
        
        loc = 5
        df = pd.read_csv(mutation_file,index_col=0)        
        plot_mutation_changes(df,gs,loc)


#         plt.tight_layout()
        plt.savefig(outdir+'/{}_{}.png'.format(bindingtype,cancertype),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
        plt.close()




    