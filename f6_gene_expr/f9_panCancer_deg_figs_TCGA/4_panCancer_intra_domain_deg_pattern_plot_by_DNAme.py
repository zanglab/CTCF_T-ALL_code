import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=17
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
#sns.set(font_scale=2)
import association_with_genes
import association_with_regions
import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
sys.path.insert(0,os.path.abspath('modules'))
import CTCF_TALL_modules_new
from scipy import stats 
from brokenaxes import brokenaxes
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'xtick.color': 'k','axes.spines.top': True,'axes.spines.bottom': True})
   
patient_label = {'PD30':'Patient2','PD9':'Patient1'}

print('\n\nNotes: up genes: [#up promoter ctrl # total promoter ctrl] [# up promoter cancer # total promoter cancer]\
    [#up domain ctrl # total domain ctrl] [# up domain cancer # total domain cancer]')        
print('down genes, same way')


def intra_domain_diff_gene_plot_nobreak(cancertypes,pos_list,neg_list,figname,ylim1,ylim2,xticklabels=None,legend=None,ylabels=None,title=None,legend_loc=None):
   
    colors = ['w','silver','w','k']
    fig,axes = plt.subplots(nrows=2,sharex=False,figsize=(5,4))
    fig.subplots_adjust(hspace=.31)
    poses = np.arange(len(pos_list))
    for pos in poses:
        print(cancertypes[pos])
        print('upgenes',pos_list[pos])
        print('dngenes',neg_list[pos],'\n')
        a1_s,a1_t = pos_list[pos][0]# const-p-s, const-p-t,
        b1_s,b1_t = pos_list[pos][1]# df-p-s, df-p-t,
        c1_s,c1_t = pos_list[pos][2]# const-d-s, const-d-t,
        d1_s,d1_t = pos_list[pos][3]# df-p-s, df-p-t,
        #m = axes[0].bar([pos,pos+0.2,pos+0.4,pos+0.6],[100*a1_s/a1_t,100*b1_s/b1_t,100*c1_s/c1_t,100*d1_s/d1_t],width=0.18,color=colors,edgecolor = "k",lw=1)
        m = axes[0].bar([pos,pos+0.4],[100*c1_s/c1_t,100*d1_s/d1_t],width=0.3,color=colors[2:],edgecolor = "k",lw=1)
#         s1,p1 = stats.fisher_exact([[a1_s,a1_t-a1_s],[b1_s,b1_t-b1_s]])
        s2,p2 = stats.fisher_exact([[c1_s,c1_t-c1_s],[d1_s,d1_t-d1_s]])
#         if p1<0.05:
#             star_mark="*"
#             if p1<0.001:
#                 star_mark="**"
#             axes[0].text(pos+0.2,100*b1_s/b1_t,star_mark,ha='center',va='center')
        if p2<0.05:
            star_mark="*"
            if p2<0.001:
                star_mark="**"
            axes[0].text(pos+0.4,100*d1_s/d1_t,star_mark,ha='center',va='center')

        a2_s,a2_t = neg_list[pos][0]
        b2_s,b2_t = neg_list[pos][1]
        c2_s,c2_t = neg_list[pos][2]
        d2_s,d2_t = neg_list[pos][3]
        #m = axes[1].bar([pos,pos+0.2,pos+0.4,pos+0.6],[100*a2_s/a2_t,100*b2_s/b2_t,100*c2_s/c2_t,100*d2_s/d2_t],width=0.18 ,color=colors,edgecolor = "k",lw=1)
        m = axes[1].bar([pos,pos+0.4],[100*c2_s/c2_t,100*d2_s/d2_t],width=0.3 ,color=colors[2:],edgecolor = "k",lw=1)
#         s1,p1 = stats.fisher_exact([[a2_s,a2_t-a2_s],[b2_s,b2_t-b2_s]])
        s2,p2 = stats.fisher_exact([[c2_s,c2_t-c2_s],[d2_s,d2_t-d2_s]])

        #if p1<0.05 and 100*b2_s/b2_t<ylim2:
#         if p1<0.05:
#             star_mark="*"
#             if p1<0.001:
#                 star_mark="**"
#             axes[1].text(pos+0.2,100*b2_s/b2_t,star_mark,ha='center',va='top')
        #if p2<0.05 and 100*d2_s/d2_t<ylim2:
        if p2<0.05:
            star_mark="*"
            if p2<0.001:
                star_mark="**"
            axes[1].text(pos+0.4,100*d2_s/d2_t,star_mark,ha='center',va='top')
    
    if legend_loc=="lower right":
        axes[1].legend([m[0],m[1]],legend,frameon=False,fontsize=14,loc="lower right",borderaxespad=-1.5,labelspacing=.2,handlelength=1) 
    else:
        axes[0].legend([m[0],m[1]],legend,frameon=False,fontsize=14,loc="upper right",borderaxespad=-1.5,labelspacing=.2,handlelength=1) 
    
    axes[0].set(xticks = poses+0.3,xticklabels = ['\n'.join(i.split('_')) for i in xticklabels],ylabel=ylabels[0])
    axes[0].set_title(title,va='baseline')
    axes[0].xaxis.set_ticks_position('bottom')
    axes[1].set(ylabel=ylabels[1])
    vals = axes[0].get_yticks()
    axes[0].set_ylim([0,ylim1])
    axes[1].set_ylim([0,ylim2])
    vals = axes[1].get_yticks()
    axes[1].set_xticks([])
    axes[0].set_yticks([0,ylim1])
    axes[1].set_yticks([0,ylim1])
    #axes[1].xaxis.tick_top()
#     sns.despine(offset=0, trim=False)
    axes[0].tick_params(axis='y',direction='out', length=3, width=.8)
    axes[0].tick_params(axis='x',direction='out', length=0, width=.8, colors='k')
    axes[1].tick_params(axis='y',direction='out', length=3, width=.8, colors='black')
    #fig.tight_layout(pad=0.1, w_pad=0.1, h_pad=.1)
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    axes[1].spines['bottom'].set_visible(False)
    axes[1].invert_yaxis()
    
    #axes[0].spines['top'].set_visible(False)
    #axes[0].spines['top'].set_visible(False)
    #axes[1].spines['top'].set_visible(False)
    #axes[1].spines['top'].set_visible(False)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.05,dpi=600,transparent=True)
    plt.close()
    


def return_intra_deg_list(df,const_df):

    pcs = ['promoter_allgenes','promoter_upgenes','promoter_dngenes'] # promoter columns
    dcs = ['domain_allgenes','domain_upgenes','domain_dngenes'] # domain columns
    
    a1,a2,a3,a4 = \
    [const_df[pcs[1]].sum(),const_df[pcs[0]].sum()],\
    [df[pcs[1]].sum(),df[pcs[0]].sum()],\
    [const_df[dcs[1]].sum(),const_df[dcs[0]].sum()],\
    [df[dcs[1]].sum(),df[dcs[0]].sum()]
    
    df[dcs[1]].sum(),df[dcs[0]].sum()
    b1,b2,b3,b4 = \
    [const_df[pcs[2]].sum(),const_df[pcs[0]].sum()],\
    [df[pcs[2]].sum(),df[pcs[0]].sum()],\
    [const_df[dcs[2]].sum(),const_df[dcs[0]].sum()],\
    [df[dcs[2]].sum(),df[dcs[0]].sum()]
    
    return [a1,a2,a3,a4],[b1,b2,b3,b4]



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

    return gained_df,lost_df,const_df


 
def main_sub(suf_name,methy_type):

    outdir = 'f4_panCancer_intra_domain_deg_pattern_figs_by_DNAme'+os.sep+suf_name
    os.makedirs(outdir,exist_ok=True)

    # intra-domain deg
    intra_domain_deg_pos_gained_list,intra_domain_deg_neg_gained_list=[],[]
    intra_domain_deg_pos_lost_list,intra_domain_deg_neg_lost_list=[],[]
    
    # ==== cancer specific CTCFs
    cancertypes=['T-ALL','AML','BRCA','CRC','LUAD','PRAD','PRAD_TissueAdded']
    cancertypes=['T-ALL','AML','BRCA','CRC','LUAD','PRAD_TissueAdded']
    cancertypes_ticks=['T-ALL','AML','BRCA','CRC','LUAD','PRAD']
    cancertypes=['T-ALL','BRCA','CRC','LUAD']
    cancertypes_ticks=['T-ALL','BRCA','CRC','LUAD']

    for cancertype in cancertypes:
        gained_df,lost_df,const_df = return_cancer_specific_binding(suf_name,cancertype)
        columns = ['differential_methy_EachSide150bp','DE_count_EachSide150bp']
        if methy_type == "methylation increased":
            lost_df = lost_df.loc[(lost_df[columns[0]]>20)&(lost_df[columns[1]]>5)]
        else:
            lost_df = lost_df.loc[(lost_df[columns[0]]<=20)|(lost_df[columns[1]]<=5)]
            
        ###
        a,b = return_intra_deg_list(gained_df,const_df)#;print(cancertype);print('gained\n',a,b)
        intra_domain_deg_pos_gained_list.append(a)
        intra_domain_deg_neg_gained_list.append(b)

        a,b = return_intra_deg_list(lost_df,const_df)#;print(cancertype);print('lost\n',a,b)
        intra_domain_deg_pos_lost_list.append(a)
        intra_domain_deg_neg_lost_list.append(b)


    figname = outdir+os.sep+'intra_domain_diff_gene_enrichment_gained_selected_ctrl_patient.pdf'
    xticklabels= cancertypes_ticks
#     legend=['promoter ctrl','promoter','domain ctrl','intra domain']
    legend=['domain ctrl','intra domain']
    ylabels=['% genes \n up-regulated','% genes \ndown-regulated']
    title='Genes near gained CTCF\n'  
#     ylim1,ylim2= 35,35
    #ylim1,ylim2 = [30,40,60],40
#     print('\n===gained-fig values===\n')
#     intra_domain_diff_gene_plot_nobreak(cancertypes,intra_domain_deg_pos_gained_list,intra_domain_deg_neg_gained_list,figname,ylim1,ylim2,xticklabels=xticklabels,legend=legend,ylabels=ylabels,title=title,legend_loc="upper right")    
       
    figname = outdir+os.sep+'intra_domain_diff_gene_enrichment_lost_selected_ctrl_patient_{}.pdf'.format('_'.join(methy_type.split(' ')))
    title='Genes near lost CTCF\n{}\n'.format(methy_type) 
    #ylim1,ylim2 = 50,[45,60,66]
    ylim1,ylim2= 65,65
    print('\n===lost-fig values===\n')
    intra_domain_diff_gene_plot_nobreak(cancertypes,intra_domain_deg_pos_lost_list,intra_domain_deg_neg_lost_list,figname,ylim1,ylim2,xticklabels=xticklabels,legend=legend,ylabels=ylabels,title=title,legend_loc="upper right")    
        

def main():

    suf_names = ['domain_GT100K_LT1M_log2FC0.585_padj1e-3','domain_GT100K_LT1M_log2FC1_padj1e-5',\
                 'domain_NoLimit_log2FC0.585_padj1e-3','domain_NoLimit_log2FC1_padj1e-5']
#     suf_names = ['domain_GT100K_LT1M_log2FC1_padj1e-5']

    for suf_name in suf_names:
        main_sub(suf_name,'methylation increased')   
        main_sub(suf_name,'methylation NOT increased')    



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

