import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=17
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.0)
#sns.set_style("whitegrid", {'axes.grid' : False})
import association_with_genes
import association_with_regions
import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
sys.path.insert(0,os.path.abspath('modules'))
import CTCF_TALL_modules_new
from scipy import stats 
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
  

def intra_domain_hicor_gene_plot(cancertypes,pos_list,plot_diff_list,figname,ylim1,ylim2,xticklabels=None,legend=None,ylabels=None,title=None):
    star_marker='*'
    colors = ['lightgreen','skyblue']
    colors = ['silver','k']
    fig,axes = plt.subplots(nrows=2,sharex=False,figsize=(6.,4))#,gridspec_kw = {'height_ratios':[2, 2]})
    fig.subplots_adjust(hspace=.35)
    poses = np.arange(len(pos_list))

    print('hicor fig: [# hicor promoter genes # total promoter genes ] [# hicor domain genes # total domain genes]')
    print('diff gene fig: [# up hicor-promoter genes # down hicor-promoter genes  #total hicor-promoter genes] [ # up hicor domain genes # down hicor domain genes #total hicor-domain genes]')

    for pos in poses:
        print(xticklabels[pos])
        print('hicor fig:',pos_list[pos])
        print('diff gene:',plot_diff_list[pos],'\n')
        
       
        a1_s,a1_t = pos_list[pos][0]# const-p-s, const-p-t,
        b1_s,b1_t = pos_list[pos][1]# df-p-s, df-p-t,
        m = axes[0].bar([pos,pos+0.4],[100*a1_s/a1_t,100*b1_s/b1_t],width=0.32,color=colors,edgecolor = "k",lw=1)

        a0_s,a0_t = pos_list[0][0]# const-p-s, const-p-t,
        b0_s,b0_t = pos_list[0][1]# df-p-s, df-p-t,
        
        s1,p1 = stats.fisher_exact([[a1_s,a1_t-a1_s],[a0_s,a0_t-a0_s]])
        s2,p2 = stats.fisher_exact([[b1_s,b1_t-b1_s],[b0_s,b0_t-b0_s]])
        if p1<0.05:
            star_marker='*'
            if p1<0.001:
                star_marker='**'
            axes[0].text(pos,100*a1_s/a1_t,star_marker,ha='center',va='center')
        if p2<0.05:
            star_marker='*'
            if p2<0.001:
                star_marker='**'
            axes[0].text(pos+0.4,100*b1_s/b1_t,star_marker,ha='center',va='center')

        a2_up,a2_dn,a2_t = plot_diff_list[pos][0]# const-p-s, const-p-t,
        b2_up,b2_dn,b2_t = plot_diff_list[pos][1]# df-p-s, df-p-t,
        m1 = axes[1].bar([pos,pos+0.4],[100*a2_up/a2_t,100*b2_up/b2_t],width=0.32,color='red')
        m2 = axes[1].bar([pos,pos+0.4],[-100*a2_dn/a2_t,-100*b2_dn/b2_t],width=0.32,color='blue')

        a0_up,a0_dn,a0_t = plot_diff_list[0][0]# const-p-s, const-p-t,
        b0_up,b0_dn,b0_t = plot_diff_list[0][1]# df-p-s, df-p-t,
        s1,p1 = stats.fisher_exact([[a2_up,a2_t-a2_up],[a0_up,a0_t-a0_up]])
        s2,p2 = stats.fisher_exact([[b2_up,b2_t-b2_up],[b0_up,b0_t-b0_up]])
        s3,p3 = stats.fisher_exact([[a2_dn,a2_t-a2_dn],[a0_dn,a0_t-a0_dn]])
        s4,p4 = stats.fisher_exact([[b2_dn,b2_t-b2_dn],[b0_dn,b0_t-b0_dn]])
        if p1<0.05:
            star_marker='*'
            if p1<0.001:
                star_marker='**'
            axes[1].text(pos,100*a2_up/a2_t,star_marker,ha='center',va='center')
        if p2<0.05:
            star_marker='*'
            if p2<0.001:
                star_marker='**'
            axes[1].text(pos+0.4,100*b2_up/b2_t,star_marker,ha='center',va='center')
        if p3<0.05:
            star_marker='*'
            if p3<0.001:
                star_marker='**'
            axes[1].text(pos,-100*a2_dn/a2_t,star_marker,ha='center',va='top')
        if p4<0.05:
            star_marker='*'
            if p4<0.001:
                star_marker='**'
            axes[1].text(pos+0.4,-100*b2_dn/b2_t,star_marker,ha='center',va='top')
    
    axes[0].legend([m[0],m[1]],legend,loc='upper left',frameon=False,fontsize=12,borderaxespad=-.1,labelspacing=.1,handletextpad=0.2,handlelength=1) 
    axes[1].legend([m1,m2],['Up','Down'] ,loc='lower left',frameon=False,fontsize=12,borderaxespad=-.1,labelspacing=.1,handletextpad=0.2,handlelength=1) 
    axes[0].set(xticks = poses+0.2,xticklabels = xticklabels,ylabel=ylabels[0])
    #axes[0].set(xticks = poses+0.2,xticklabels = ['\n'.join(i.split('_')) for i in xticklabels],title=title,ylabel=ylabels[0])
    axes[1].set(ylabel=ylabels[1])
    axes[1].set_ylim([-100,100])
    axes[0].set_ylim([0,45])
    axes[0].set_title(title,va='bottom',fontsize=17)
    axes[1].axhline(y=0,color='k',linestyle='-',linewidth=1)
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    axes[1].spines['bottom'].set_visible(False)
    axes[1].spines['left'].set_visible(True)
    axes[1].set_xticks([])
    axes[0].tick_params(axis='y',direction='out', length=3, width=.8, colors='black')
#     axes[0].tick_params(axis='x',direction='out', length=0, width=.8, colors='black')
    axes[1].tick_params(axis='y',direction='out', length=3, width=.8, colors='black')
    
    #fig.tight_layout(pad=0.1, w_pad=0.1, h_pad=.1)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,transparent=True)
    plt.close()

def return_intra_hicor_list(df,const_df):
        
    pcs = ['promoter_allgenes','promoter_allgenes_hicor','promoter_upgenes_hicor','promoter_dngenes_hicor']  
    dcs = ['domain_allgenes','domain_allgenes_hicor','domain_upgenes_hicor','domain_dngenes_hicor']
 
    a1,a2,a3,a4 = \
    [const_df[pcs[1]].sum(),const_df[pcs[0]].sum()],\
    [df[pcs[1]].sum(),df[pcs[0]].sum()],\
    [const_df[dcs[1]].sum(),const_df[dcs[0]].sum()],\
    [df[dcs[1]].sum(),df[dcs[0]].sum()]
    #print(a1,a2,a3,a4)
    return a1,a2,a3,a4

def return_intra_hicorDiff_list(df,const_df):
        
    pcs = ['promoter_allgenes','promoter_allgenes_hicor','promoter_upgenes_hicor','promoter_dngenes_hicor']  
    dcs = ['domain_allgenes','domain_allgenes_hicor','domain_upgenes_hicor','domain_dngenes_hicor']
 
    a1,a2,a3,a4 = \
    [const_df[pcs[2]].sum(),const_df[pcs[3]].sum(),const_df[pcs[1]].sum()],\
    [df[pcs[2]].sum(),df[pcs[3]].sum(),df[pcs[1]].sum()],\
    [const_df[dcs[2]].sum(),const_df[dcs[3]].sum(),const_df[dcs[1]].sum()],\
    [df[dcs[2]].sum(),df[dcs[3]].sum(),df[dcs[1]].sum()]
    #print(a1,a2,a3,a4)
    return a1,a2,a3,a4

def return_cancer_specific_binding(suf_name,cancertype):
    # == cancer specific gained/lost bindings ====
    pardir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f7_panCancer_deg_figs/f3_FeatureCombinaiton_CancerSpecificBinding_local'
#     cancertypes=['T-ALL','AML','BRCA','CRC','LUAD','PRAD','PRAD_TissueAdded']
    gained_file = pardir+os.sep+'f3_cancer_specific_binding_{}'.format(suf_name)+os.sep+'{}_gained.csv'.format(cancertype)
    lost_file = pardir+os.sep+'f3_cancer_specific_binding_{}'.format(suf_name)+os.sep+'{}_lost.csv'.format(cancertype)
    const_file = pardir+os.sep+'f3_cancer_specific_binding_{}'.format(suf_name)+os.sep+'{}_const.csv'.format(cancertype)
    with open(gained_file) as gained_inf, open(lost_file) as lost_inf, open(const_file) as const_inf:
        gained_df = pd.read_csv(gained_inf,index_col=0)
        lost_df = pd.read_csv(lost_inf,index_col=0)
        const_df = pd.read_csv(const_inf,index_col=0)
    gained_df_new,lost_df_new = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype)
    gained_df = gained_df.loc[gained_df_new.index]
    lost_df = lost_df.loc[lost_df_new.index]
    print(cancertype,'gained',gained_df.shape)
    print(cancertype,'lost',lost_df.shape)
    print(cancertype,'const',const_df.shape)
    return gained_df,lost_df,const_df

  
def main_sub(suf_name):

    outdir = 'f6_panCancer_intra_domain_hicor_figs'+os.sep+suf_name
    os.makedirs(outdir,exist_ok=True)

    # intra-domain deg
    intra_domain_hicor_gained_list=[]
    intra_domain_hicor_lost_list=[]

    intra_domain_hicorDiff_gained_list=[]
    intra_domain_hicorDiff_lost_list=[]
    
    cancertypes=['T-ALL','AML','BRCA','CRC','LUAD','PRAD_TissueAdded']
    cancertype_sticklabels=['T-ALL','AML','BRCA','CRC','LUAD','PRAD','PRAD']
    
    for cancertype in cancertypes:
        gained_df,lost_df,const_df = return_cancer_specific_binding(suf_name,cancertype)
        
        ctrl_p,treat_p,ctrl_d,treat_d = return_intra_hicor_list(gained_df,const_df)
        intra_domain_hicor_gained_list.append([treat_p,treat_d])
        ctrl_diffp,treat_diffp,ctrl_diffd,treat_diffd = return_intra_hicorDiff_list(gained_df,const_df)
        intra_domain_hicorDiff_gained_list.append([treat_diffp,treat_diffd])
        
        ctrl2_p,treat_p,ctrl2_d,treat_d = return_intra_hicor_list(lost_df,const_df)
        intra_domain_hicor_lost_list.append([treat_p,treat_d])
        ctrl2_diffp,treat_diffp,ctrl2_diffd,treat_diffd = return_intra_hicorDiff_list(lost_df,const_df)
        intra_domain_hicorDiff_lost_list.append([treat_diffp,treat_diffd])
        
        
    intra_domain_hicor_gained_list.insert(0,[ctrl_p,ctrl_d])
    intra_domain_hicorDiff_gained_list.insert(0,[ctrl_diffp,ctrl_diffd])

    figname = outdir+os.sep+'intra_domain_hicor_genes_gained.pdf'
    xticklabels= np.append('Control',cancertype_sticklabels)
    legend=['Promoter','Intra-domain']
    ylabels=['% highly correlated \nCTCF-gene pairs','% differentially \nexpressed genes']
    title='Genes near gained CTCF'  
    ylim1,ylim2= 25.5,25 
    print('\n===gained-fig values===\n')
    intra_domain_hicor_gene_plot(cancertypes,intra_domain_hicor_gained_list,intra_domain_hicorDiff_gained_list,figname,ylim1,ylim2,xticklabels=xticklabels,legend=legend,ylabels=ylabels,title=title)
    
    
    intra_domain_hicor_lost_list.insert(0,[ctrl2_p,ctrl2_d])
    intra_domain_hicorDiff_lost_list.insert(0,[ctrl2_diffp,ctrl2_diffd])
    figname = outdir+os.sep+'intra_domain_hicor_genes_lost.pdf'
    title='Genes near lost CTCF' 
    print('\n===lost-fig values===\n')
    intra_domain_hicor_gene_plot(cancertypes,intra_domain_hicor_lost_list,intra_domain_hicorDiff_lost_list,figname,ylim1,ylim2,xticklabels=xticklabels,legend=legend,ylabels=ylabels,title=title)
    
     
def main():

    suf_names = ['domain_GT100K_LT1M_log2FC0.585_padj1e-3','domain_GT100K_LT1M_log2FC1_padj1e-5',\
                 'domain_NoLimit_log2FC0.585_padj1e-3','domain_NoLimit_log2FC1_padj1e-5']
#     suf_names = ['domain_GT100K_LT1M_log2FC1_padj1e-5']

    for suf_name in suf_names:
        main_sub(suf_name)    
   

    
    
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

