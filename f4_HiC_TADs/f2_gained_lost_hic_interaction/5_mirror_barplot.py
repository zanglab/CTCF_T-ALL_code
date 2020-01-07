import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=18
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
# import myplot
import CTCF_TALL_modules_new
matplotlib.rcParams["font.sans-serif"] = ["Arial"]

chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']

cancer_type_hic_filename = {'T-ALL':['Jurkat','A6010'],\
                                'T-ALL-P1':['PD9','A6010'],\
                                'T-ALL-P2':['PD31','A6010'],\
                                'CRC':['HCT116','trans_colon1'],\
                                'PRAD':['LNCaP','PrEC'],\
                                'PRAD_TissueAdded':['LNCaP','PrEC']}

cancertype_title = {'T-ALL':'T-ALL','T-ALL-P1':'T-ALL patient1','T-ALL-P2':'T-ALL patient2',\
                    'CRC':'CRC','PRAD':'PRAD','PRAD_TissueAdded':'PRAD'}

cancertype_match_title = {'T-ALL':'T-ALL vs. T cell',\
                          'T-ALL-P1':'T-ALL patient1 vs. T cell',\
                          'T-ALL-P2':'T-ALL patient2 vs. T cell',\
                          'CRC':'CRC vs. normal colon tissue',\
                          'PRAD':'PRAD vs. prostate tissue',\
                          'PRAD_TissueAdded':'PRAD vs. prostate tissue'}

pvalue_type,stats_type = 'pvalue','stats'
#pvalue_type,stats_type = 'log_pvalue','log_stats'


#########
#########

def mirror_bar_plot_with_fisher_test(cancertype,plot_list,plot_list2,colors,figname,xticklabels=None,legend=None,ylabels=None,title=None):

    fig,axes = plt.subplots(nrows=2,sharex=True,figsize=(6,4))
    poses = np.arange(len(plot_list))
    for pos in poses:
        a1_s,a1_t = plot_list[pos][0]
        b1_s,b1_t = plot_list[pos][1]
        c1_s,c1_t = plot_list[pos][2]
        m = axes[0].bar([pos,pos+0.23,pos+0.46],[a1_s/a1_t,b1_s/b1_t,c1_s/c1_t],width=0.22,color=colors)
        s1,p1 = stats.fisher_exact([[b1_s,b1_t-b1_s],[a1_s,a1_t-a1_s]])
        s2,p2 = stats.fisher_exact([[c1_s,c1_t-c1_s],[a1_s,a1_t-a1_s]])
        p_mark1 = '*'
        p_mark2 = '*'
        if p1<0.05 :
            if p1<0.001:
                p_mark1='**'
            axes[0].text(pos+0.23,b1_s/b1_t,p_mark1,ha='center',va='center')
        if p2<0.05 :
            if p2<0.001:
                p_mark2='**'
            axes[0].text(pos+0.46,c1_s/c1_t,p_mark2,ha='center',va='center')

        a2_s,a2_t = plot_list2[pos][0]
        b2_s,b2_t = plot_list2[pos][1]
        c2_s,c2_t = plot_list2[pos][2]
        m = axes[1].bar([pos,pos+0.23,pos+0.46],[a2_s/a2_t,b2_s/b2_t,c2_s/c2_t],width=0.22,color=colors)
        s1,p1 = stats.fisher_exact([[b2_s,b2_t-b2_s],[a2_s,a2_t-a2_s]])
        s2,p2 = stats.fisher_exact([[c2_s,c2_t-c2_s],[a2_s,a2_t-a2_s]])
        p_mark1 = '*'
        p_mark2 = '*'
        if p1<0.05 :
            if p1<0.001:
                p_mark1='**'
            axes[1].text(pos+0.23,b2_s/b2_t,p_mark1,ha='center',va='top')
        if p2<0.05  :
            if p2<0.001:
                p_mark2='**'
            axes[1].text(pos+0.46,c2_s/c2_t,p_mark2,ha='center',va='top')
    
    fig.legend([m[0],m[1],m[2]],legend,bbox_to_anchor=[1.28,.95],frameon=False,fontsize=14,handlelength=1) 
    axes[0].set(xticks = poses+0.1,xticklabels = xticklabels,ylabel=ylabels[0])
    axes[0].set_title(title+'\n',va='center')
    axes[1].set(ylabel=ylabels[1])
#     axes[0].set_ylim([0,0.60])
#     axes[1].set_ylim([0,0.60])
    #axes[0].set_xticklabels(xticklabels,ha='center',rotation=45)
    vals = axes[0].get_yticks()
    axes[0].set_yticklabels(['{:,.0%}'.format(x) for x in vals])
    vals = axes[1].get_yticks()
    axes[1].set_yticklabels(['{:,.0%}'.format(x) for x in vals])
    axes[1].xaxis.tick_top()
    axes[1].invert_yaxis()
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['bottom'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    fig.tight_layout(pad=0.1, w_pad=0.1, h_pad=.1)
    plt.savefig(figname,bbox_inches='tight',pad_inches=.05,transparent=True,dpi=600)
    plt.close()


def return_sig_nums(cancertype,df,pvalue_thre=0.05):

    total = len(df.index)
    pos_num = len(df.loc[(df[stats_type]>0)&(df[pvalue_type]<pvalue_thre)].index)
    minus_num = len(df.loc[(df[stats_type]<0)&(df[pvalue_type]<pvalue_thre)].index)
    
    df = df.loc[df[pvalue_type]<0.05]
    cancer_mean = '{}_mean'.format(cancer_type_hic_filename[cancertype][0])
    normal_mean = '{}_mean'.format(cancer_type_hic_filename[cancertype][1])
    values = np.log10(df[cancer_mean]/df[normal_mean])
    values = values.replace([np.inf,-np.inf],np.nan).dropna(0)
    
    return total,pos_num,minus_num,values



def read_ttest_df(indir,cancertype,normalization,resolution,viewregion):

    gained_df = pd.read_csv(indir+os.sep+'{}_{}_{}_res-{}_viewregion-{}_paired_ttest.csv'.format(cancertype,'gained',normalization,resolution,viewregion),index_col=0)
    gained_df = gained_df.fillna(0)
    lost_df = pd.read_csv(indir+os.sep+'{}_{}_{}_res-{}_viewregion-{}_paired_ttest.csv'.format(cancertype,'lost',normalization,resolution,viewregion),index_col=0)
    lost_df = lost_df.fillna(0)
    const_df = pd.read_csv(indir+os.sep+'{}_{}_{}_res-{}_viewregion-{}_paired_ttest.csv'.format(cancertype,'ctrl',normalization,resolution,viewregion),index_col=0)
    const_df = const_df.fillna(0)
    
    if cancertype in ['T-ALL-P1','T-ALL-P2']:
        gained_df_new,lost_df_new = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    else:
        gained_df_new,lost_df_new = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype)
    
    # ==== in case you may need to change the control bindings/ patient verified bindings
#     gained_df2,lost_df2,const_df2 = CTCF_TALL_modules.return_specific_binding_df(cancertype)
#     const_df = const_df.loc[const_df2.index]#;print(const_df.shape);exit()
    return gained_df,lost_df,const_df


def main_cancertype(cancertype):
    
    indir = 'f3_compr_ttest_csv'
    outdir = 'f5_mirror_bar_figs'
    os.makedirs(outdir,exist_ok=True)
    viewregions = [20000,50000,100000,200000,500000,1000000]
    #viewregions = [20000,50000,100000,200000,500000,1000000]
    
    for normalization in  ['raw']:
        for resolution in [5000]:
            pos_list,minus_list = [],[]
            xticklabels = []
            for viewregion in viewregions:
                if viewregion > resolution and viewregion%resolution==0:
                    xticklabels.append('{}k'.format(int(viewregion/1000)) )
                    gained_df,lost_df,const_df = read_ttest_df(indir,cancertype,normalization,resolution,viewregion)        
                    
                    gained_total,gained_pos,gained_minus,gained_values= return_sig_nums(cancertype,gained_df);#print(gained_total,gained_pos,gained_minus)
                    lost_total,lost_pos,lost_minus,lost_values = return_sig_nums(cancertype,lost_df);#print(lost_total,lost_pos,lost_minus)
                    const_total,const_pos,const_minus,const_values = return_sig_nums(cancertype,const_df) ;#print(const_total,const_pos,const_minus)        

                    pos_list.append([[const_pos,const_total],[lost_pos,lost_total],[gained_pos,gained_total]])
                    minus_list.append([[const_minus,const_total],[lost_minus,lost_total],[gained_minus,gained_total]])
                
            colors = ['grey','skyblue','salmon']
            figname = outdir+os.sep+'{}_{}_res{}_mirror_bar.pdf'.format(cancertype,normalization,resolution)
            legend = ['Constitutive','Lost','Gained']
            ylabels = ['Increased Hi-C\ninteractions (%)','Decreased Hi-C\ninteractions (%)']
            title = '{}'.format(cancertype_match_title[cancertype])
            mirror_bar_plot_with_fisher_test(cancertype,pos_list,minus_list,legend = legend,colors=colors,figname=figname,xticklabels = xticklabels,ylabels=ylabels,title=title);
                

def main():

    cancertypes=['T-ALL','T-ALL-P1','T-ALL-P2','CRC']#,'PRAD','PRAD_TissueAdded']
    for cancertype in cancertypes:
        main_cancertype(cancertype)#;exit()
               
                
                


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--viewregion', action = 'store', type = int,dest = 'viewregion', help = 'input file of', metavar = '<int>')
    parser.add_argument('-n', '--normalization', action = 'store', type = str,dest = 'normalization', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-c', '--chrom', action = 'store', type = str,dest = 'chrom', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
