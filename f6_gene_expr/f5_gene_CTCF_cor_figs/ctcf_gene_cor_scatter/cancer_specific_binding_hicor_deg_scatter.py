'''
This file is used to get nearby geneInfo for each CTCF binding
'''

import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from GenomeData import *
from scipy import stats
import association_with_regions
import re,bisect
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=15
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import CTCF_TALL_modules_new
import json

# ==== intra-domain hicor gene
json_file_dir= '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/f3_CTCF_domain_allGene_hicorGene'
#     with open(json_file_dir+os.sep+'CTCF_domain_NoLimit_allGene_hicorGene.json', 'r') as json_file:
with open(json_file_dir+os.sep+'CTCF_domain_GT100K_LT1M_allGene_hicorGene.json', 'r') as json_file:
    ctcf_domain_gene_dic = json.load(json_file)        

# ==== CTCF binding
binding_signal_sqrt_QN_file = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f3_CTCF_chip_RPKM/f1_chip_RPKM_on_union_bindings/RPKM_on_all_union_bindings_RPKM_CellType_combined_sqrt_QN.csv"
with open(binding_signal_sqrt_QN_file) as binding_signal_sqrt_QN_inf:
    binding_signal = pd.read_csv(binding_signal_sqrt_QN_inf,index_col = 0)
    
# ==== gene expression
gene_expr_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f2_gene_expr_TPM/f5_gene_expr_combined/geneExpr_replicates_combined_TPM_sqrt.csv"
with open(gene_expr_file) as gene_expr_inf:
    gene_expr = pd.read_csv(gene_expr_inf,index_col = 0)    



def scatter_plot_ctcf_gene(ctcf_id,gene_id,figname):    
    #print(sorted(binding_signal.columns));exit(0)
    samples = [ii for ii in binding_signal.columns]
    binding_val = binding_signal.loc[ctcf_id][samples]
    expr_val = gene_expr.loc[gene_id][samples]
    # calculate the R score
    slope, intercept, r_value, p_value, std_err = stats.linregress(binding_val, expr_val)       
    # plot the scatter plot        
    plt.figure(figsize=(3,3))
    for sample in samples:
        n = plt.scatter(binding_val[sample],expr_val[sample],c = 'navy')
    # plot the fitted line
    x = np.sort(binding_val)
    plt.plot(x,x*slope+intercept,c = 'grey',ls='--',lw=.6)
    plt.text(.65,.02,'$R^2$ = {:.2f}'.format(r_value**2),fontsize=12,transform=plt.axes().transAxes)
    #plt.legend([n,c],['Tissue','Cell Line'],frameon=True,fontsize=10,loc=0) 
    plt.xlabel('CTCF binding (sqrt RPKM)')
    plt.ylabel('gene expression (sqrt TPM)')
#     plt.title('binding{} ~ {}'.format(ctcf_id,gene_id))
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = 0.1,transparent=True)
    plt.close()



def scatter_each_cancertype(cancertype,outdir):
    
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype)
    #print(gained_df.columns)
    for ctcf_id in gained_df.index:
        hicor_gene = ctcf_domain_gene_dic[str(ctcf_id)]['intra_domain_genes_hicor']
        if len(hicor_gene)>0:
            for gene_id in hicor_gene:
                figname = '{}/gained_hicor_gene_{}_{}_{}.pdf'.format(outdir,cancertype,ctcf_id,gene_id)
                scatter_plot_ctcf_gene(ctcf_id,gene_id,figname)

    for ctcf_id in lost_df.index:
        hicor_gene = ctcf_domain_gene_dic[str(ctcf_id)]['intra_domain_genes_hicor']
        if len(hicor_gene)>0:
            for gene_id in hicor_gene:
                figname = '{}/lost_hicor_gene_{}_{}_{}.pdf'.format(outdir,cancertype,ctcf_id,gene_id)
                scatter_plot_ctcf_gene(ctcf_id,gene_id,figname)


def main():

    outdir = 'cancer_specific_binding_hicor_deg_scatter'
    os.makedirs(outdir,exist_ok=True)

    cancertypes=['T-ALL','AML','BRCA','CRC','LUAD','PRAD_TissueAdded']
#     cancertypes=['T-ALL']
    for cancertype in cancertypes:
        scatter_each_cancertype(cancertype,outdir)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--geneid', action = 'store', type = str,dest = 'geneid', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
