import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import scipy
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.4)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")

import CTCF_TALL_modules_new
matplotlib.rcParams["font.sans-serif"] = ["Arial"]



# def return_ctcf_intra_domain_genes(binding_type,colname):
#     # intra domain genes of CTCF bindings
#     ctcf_domain_gene_dic = CTCF_TALL_modules_new.return_union_domain_gene()
#     gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
#     
#     tall_gained_intra_domain_genes = np.array([])
#     for gained_id in gained_df.index:
#         domain_genes = ctcf_domain_gene_dic[str(gained_id)]['intra_domain_genes']
#         tall_gained_intra_domain_genes = np.append(tall_gained_intra_domain_genes,domain_genes)
#     return tall_gained_intra_domain_genes


def return_essential_gene_df(df,cancertype):

    if cancertype =='BRCA':
        columns = ['T47D_Day21_vs_T47D_Day0|neg|score','T47D_Day21_vs_T47D_Day0|neg|p-value','T47D_Day21_vs_T47D_Day0|neg|fdr']
    else:
        columns = ['LNCaP_Day21_vs_LNCaP_Day0|neg|score','LNCaP_Day21_vs_LNCaP_Day0|neg|p-value','LNCaP_Day21_vs_LNCaP_Day0|neg|fdr']
                    
    # ==== essential genes as genes pass fdr or p threshold
    fdr_thre = 0.25
    p_thre=0.05
#     df_essential = df[df[columns[1]]<fdr_thre]
#     df_essential = df[df[columns[0]]<p_thre]
    df_essential = df.sort_values(by=[columns[0]]).iloc[:int(df.shape[0]*.1)]
    return df_essential


def plot_compr_bar(bars,pvalue,colors,xticklabels,basename):

    # ====  plot
    plt.figure(figsize=(1.5,2.8))
    poses = np.arange(len(xticklabels))
    plt.bar(poses,bars,color=colors,width=0.65,lw=0)
                    
    if pvalue<0.05:
        star_mark="*"
        if pvalue<0.001:
            star_mark="**"
        plt.text(poses[-1],bars[-1],star_mark,ha='center',fontsize=17)

#     sns.despine(offset=0, trim=False)
    plt.ylabel('% essential genes',fontsize=18)
    plt.ylim([0,22])
    plt.xlim([-0.6,1.6])
    plt.axes().set_xticks(poses)
    plt.axes().set_xticklabels(xticklabels,rotation=30)
    plt.savefig(basename+'_compr.pdf',transparent=True,bbox_inches='tight',pad_inches=0.01,dpi=600)
    plt.close()
                   



def main():


    outdir = 'f4_hicor_AllGene_essentiality_figs_layout'
    os.makedirs(outdir,exist_ok=True)
    
    crispr_results_file = "../data/supp_tables/pnas.1908155116.sd01.xlsx"
    crispr_results = pd.read_excel(crispr_results_file,index_col=0)
    
    ctcf_domain_gene_dic = CTCF_TALL_modules_new.return_union_domain_gene()
    cancertypes=['BRCA','PRAD','PRAD_TissueAdded']
    
    for cancertype in cancertypes:
        basename = outdir+os.sep+'{}_hicor_gene'.format(cancertype)     
        
        # ==== collect intra-domain hi-cor genes
        gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype)
        gained_intra_domain_genes = np.array([])
        for gained_id in gained_df.index:
#             domain_genes = ctcf_domain_gene_dic[str(gained_id)]['intra_domain_genes']
            domain_genes = ctcf_domain_gene_dic[str(gained_id)]['intra_domain_genes_hicor']
            gained_intra_domain_genes = np.append(gained_intra_domain_genes,domain_genes)
        
        # ==== essential genes in cancer
        crispr_essential =  return_essential_gene_df(crispr_results,cancertype)   
        # ==== intra-domain hi-cor gene
        crispr_hicor = crispr_results.loc[set(gained_intra_domain_genes)].dropna()
        # ==== essential intra-domain hi-cor gene
#         crispr_hicor_essential = return_essential_gene_df(crispr_hicor,cancertype)     
        crispr_hicor_essential = crispr_essential.loc[set(gained_intra_domain_genes)].dropna() 
        crispr_hicor.to_csv(basename+'.csv') 
        crispr_hicor_essential.to_csv(basename+'_essential.csv')     
        
        colors = ['lightgrey','k']
        xticklabels = ['All genes','Intra-domain\nhi-cor genes']
        c1,c2 = crispr_essential.shape[0],crispr_results.shape[0]
        t1,t2 = crispr_hicor_essential.shape[0],crispr_hicor.shape[0]
        print(basename,'\n',c1,c2,t1,t2)
        s,p = stats.fisher_exact([[t1,t2-t1],[c1,c2-c1]]);print(s,p)  
        bars = [100*c1/c2,100*t1/t2]                  
        plot_compr_bar(bars,p,colors,xticklabels,basename)





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

