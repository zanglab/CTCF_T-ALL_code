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
import CTCF_TALL_modules
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'w'})

#plus = re.compile('\+')
#minus = re.compile('\-')

# cancer_type_matched = {'T-ALL':['CD4','JURKAT','CUTLL1','PD9','PD31','PD40','PTBG'],'Breast_cancer':['Breast_normal','Breast_cancer'],'Lung_cancer':['Lung_normal','Lung_cancer',],'Colon_cancer':['Colon_normal','Colon_cancer',]}
# cancer_type_matched_title = {'T-ALL':['CD4','Jurkat','Cutll1','PD9','PD31','PD40','PTBG'],'Breast_cancer':['Breast \ntissue','Breast \ncancer',],'Lung_cancer':['Lung \ntissue','Lung \ncancer',],'Colon_cancer':['Colon \ntissue','Colon \ncancer',]}

cancer_type_matched = {'T-ALL':['CD4','JURKAT'],'Breast_cancer':['Breast_normal','Breast_cancer'],'Lung_cancer':['Lung_normal','Lung_cancer',],'Colon_cancer':['Colon_normal','Colon_cancer',]}
cancer_type_matched_title = {'T-ALL':['T-cell','T-ALL'],'Breast_cancer':['Breast \nTissue','BRCA',],'Lung_cancer':['Lung \nTissue','LUAD',],'Colon_cancer':['Colon \nTissue','CRC',]}

cancer_type_matched_methylation = {'Lung_cancer':['lung','A549',],'Colon_cancer':['sigmoid_colon','HCT116',],'Breast_cancer':['breast','MCF7',],\
                        'T-ALL':['CD4','Jurkat'],'PD9':['CD4','PD9',],'PD31':['CD4','PD31',],'PD40':['CD4','PD40',],'PTBG':['CD4','PTBG',]}


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
        


def cancertype_heatmap_color(cancertype):
    if cancertype =='T-ALL':
        return 'orangered'
    if cancertype =='Breast_cancer':
        return 'goldenrod'
    if cancertype =='Colon_cancer':
        return 'darkcyan'
    if cancertype =='Lung_cancer':
        return 'royalblue'


def window_smooth(df,half_window):
    smooth_df = pd.DataFrame(index=df.index,columns=df.columns)
    for col in np.arange(len(df.columns)):
        window_left = max(col-half_window,0)
        window_right = min(col+half_window,len(df.columns)-1)
        smooth_df.iloc[:,col] = df.iloc[:,window_left:window_right+1].mean(axis=1)
    return smooth_df   
    print(df,smooth_df);exit()
    return df

def return_methylation_state_df(cancertype,samplename,flag):
    # return cancer/normal methylation file, for gained or lost bindings 
    methylation_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/fu_all_feature_analysis_panCancer/f9_combine_heatmap_CTCF_DNAme_HM/DNAme_heatmap_reindexed_csv'
    methylation_file = methylation_dir+os.sep+'{}_{}_{}_DNA_methylation.csv'.format(cancertype,samplename,flag)
    basename = os.path.basename(methylation_file)
    with open(methylation_file) as methylation_inf:
        methylation_df = pd.read_csv(methylation_inf,index_col=0)
    return methylation_df,basename

def return_CTCF_state_df(cancertype,cancer_subname,flag):
    # return the df for heatmap
    binding_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/fu_all_feature_analysis_panCancer/f9_combine_heatmap_CTCF_DNAme_HM/ctcf_heatmap_reindexed_csv'
    binding_file = binding_dir+os.sep+'{}_{}.csv'.format(flag,cancer_subname)
    with open(binding_file) as binding_inf:
        binding_df = pd.read_csv(binding_inf,index_col=0)
    # methylation data for rerank of index
    
    gained_df,lost_df,const_df = CTCF_TALL_modules.return_specific_binding_df(cancertype)
    columns = ['differential_methy_EachSide150bp_cov5','DE_count_EachSide150bp_cov5']
    if cancertype=="T-ALL":
        columns = ['differential_methy_EachSide150bp_Jurkat_vs_A601x','DE_count_EachSide150bp_Jurkat_vs_A601x']
    if flag=='gained':
        rank_index = gained_df.sort_values(by=columns[0],ascending=False).index
    elif flag=='lost':
        rank_index = lost_df.sort_values(by=columns[0],ascending=False).index
    
    binding_df = binding_df.loc[rank_index]
    binding_df = window_cumulative(binding_df)
    binding_df = signal_centered(binding_df)

    return binding_df


def sub_heatmap_plot(df,gs,title,loc,total,satuation,vmin,vmax,flag,cancertype):
    # plot each heatmap panel
    all_values = [i for col in df.columns for i in df[col].values]
    df = df.clip(upper=np.percentile(all_values,satuation))
    ax = plt.subplot(gs[0,loc])
    pal = sns.light_palette(cancertype_heatmap_color(cancertype),as_cmap=True)
    pal = sns.light_palette('red',as_cmap=True)
    #pal=LinearSegmentedColormap.from_list('rg',["w", "r"]) 
    cbarvmin=0
    if loc==0:   
        g=sns.heatmap(df,ax=ax,yticklabels=False,xticklabels=True,cbar=True,cmap=pal,vmin=cbarvmin,vmax = vmax,cbar_kws={"shrink": 0.5})
        ax.set_ylabel('CTCF ChIP-seq'.format(flag.capitalize()),fontsize=18)
        xp = g.get_xticks()#;print(xp)
        ax.set_xticks([xp[0],xp[-1]])
        if cancertype=='Lung_cancer':
            ax.set_xticklabels(['-1kb','1kb'],rotation=30,fontsize=16)
        else:
            ax.set_xticklabels(['',''])
        ax.tick_params(axis='x',direction='out', length=0, width=1, colors='black')    
        ax.set_title('{}'.format(title),fontsize=18)    
        cbar = g.collections[0].colorbar    
        cbar.set_clim(vmax*.15,vmax)
        cbar.remove()
        
    elif loc==total-1:
        if 0:
            g=sns.heatmap(df,ax=ax,yticklabels=False,xticklabels=False,cbar=False,cmap=pal,vmin=cbarvmin,vmax = vmax,cbar_kws={"shrink": 0.5})
        else:
            g=sns.heatmap(df,ax=ax,yticklabels=False,xticklabels=False,cbar=True,cmap=pal,vmin=cbarvmin,vmax = vmax,cbar_kws={"shrink": 0.5})
            cbar = g.collections[0].colorbar
            cbar.set_ticks([cbarvmin,vmax])
            cbar.set_ticklabels([vmin,vmax])
            #pos = cbar.ax.get_position() 
            #pos.x0+=0.1  
            cbar.ax.set_position([.9,0.35,.8,.5])         
            cbar.set_clim(vmax*.15,vmax)
            cbar.remove()
        #ax.text(200,-1*df.shape[0]*0.07,'{}'.format(df.shape[0]),fontsize=13,ha='right')
        ax.text(225,0,'{}'.format(df.shape[0]),fontsize=18,ha='left')
        #ax.text(-10,df.shape[0]*0.95,'0',fontsize=11,ha='right')
        ax.set_title('{}'.format(title),fontsize=18)
        ax.set_ylabel('')
        
     
   

def combine_plot(cancertypes,flag,outdir,satuation,vmin,vmax):
    # for each cancer type, return cancer/normal df for heatmap
    total_len = 2

    fig = plt.figure(figsize = (3.2,2.7))
    # to make sure the last figure have more space for color bar
    width_ratio = [1,1]
 #    if cancertypes[0]=='Lung_cancer':
#         width_ratio = [1,1.22]
        
    gs = gridspec.GridSpec(1,total_len,width_ratios=width_ratio,wspace=0.0) 
    i=0
    for cancertype in cancertypes:
        for loc_index in np.arange(len(cancer_type_matched[cancertype])):
            cancer_subname = cancer_type_matched[cancertype][loc_index]
            cancer_title = cancer_type_matched_title[cancertype][loc_index]
            loc_df = return_CTCF_state_df(cancertype,cancer_subname,flag)
            #loc_df.to_csv(outdir+os.sep+'{}_CTCF_{}.matrix.csv'.format(flag,cancertypes[0]))
            sub_heatmap_plot(loc_df,gs,cancer_title,i,total_len,satuation,vmin,vmax,flag,cancertype)
            i+=1 
    plt.savefig(outdir+os.sep+'{}_CTCF_{}.png'.format(flag,cancertypes[0]),bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close();




def main():

    outdir = 'f7_CTCF_heatmap_sep_cancer_rerank_by_py6'
    os.makedirs(outdir,exist_ok=True)
    #cancertypes = ['T-ALL','Breast_cancer','Colon_cancer','Lung_cancer']
    for cancertypes in [ ['T-ALL'], ['Breast_cancer'],['Colon_cancer'],['Lung_cancer']]:
        vmin=0
        for satuation in [98]:
            for vmax in [120]:
                pass
                #combine_plot(cancertypes,'gained',outdir,satuation,vmin,vmax)#;exit()
                #combine_plot(cancertypes,'lost',outdir,satuation,vmin,vmax)
            

    vmax = 120
    fig,ax = plt.subplots(figsize=(1.,.2))
    norm = matplotlib.colors.Normalize(vmin=0, vmax=vmax)
    pal = sns.light_palette('red',as_cmap=True)
    cb = matplotlib.colorbar.ColorbarBase(ax,cmap=pal,norm=norm,orientation='horizontal')
#     cb.set_label('change level',rotation=90,labelpad=-110)
    cb.set_ticks([0,vmax])
    cb.set_clim(vmax*.15,vmax)
    ax.tick_params(axis='x',direction='out', length=0, width=1, colors='black')    
    plt.savefig(outdir+os.sep+'colorbar.png',bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()

    fig,ax = plt.subplots(figsize=(.25,1))
    norm = matplotlib.colors.Normalize(vmin=0, vmax=vmax)
    pal = sns.light_palette('red',as_cmap=True)
    cb = matplotlib.colorbar.ColorbarBase(ax,cmap=pal,norm=norm,orientation='vertical')
#     cb.set_label('change level',rotation=90,labelpad=-110)
    cb.set_ticks([0,vmax])
    cb.set_clim(vmax*.15,vmax)
    plt.savefig(outdir+os.sep+'colorbar2.png',bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()
# 
# 
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

