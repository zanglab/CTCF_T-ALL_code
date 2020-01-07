import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
from scipy import stats

import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')


def main(gene_id):
    
#     gene_id = 'MYC'
    binding_signal_sqrt_QN_file = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f3_CTCF_chip_RPKM/f1_chip_RPKM_on_union_bindings/RPKM_on_all_union_bindings_RPKM_CellType_combined_sqrt_QN.csv"
    with open(binding_signal_sqrt_QN_file) as binding_signal_sqrt_QN_inf:
        binding_signal = pd.read_csv(binding_signal_sqrt_QN_inf,index_col = 0)
    
    gene_expr_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f2_gene_expr_TPM/f5_gene_expr_combined/geneExpr_replicates_combined_TPM_sqrt.csv"
    with open(gene_expr_file) as gene_expr_inf:
        gene_expr = pd.read_csv(gene_expr_inf,index_col = 0)    
    
    samples = [ii for ii in binding_signal.columns]
    expr_val = gene_expr.loc[gene_id][samples]
    
    outdir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/f1_gene_CTCF_cor/f7_cor_CTCF_binding_geneExpr/CTCF_gene_cor_csv'
    os.makedirs(outdir,exist_ok=True)
    
    binding_expr_cor = open(outdir+os.sep+'{}_cor.csv'.format(gene_id),'w')
    #binding_expr_cor.write('{}\t{}\n'.format('rvalue','pvalue'))
    for binding_id in binding_signal.index:
        binding_val = binding_signal.loc[binding_id][samples]
        slope, intercept, r_value, p_value, std_err = stats.linregress(binding_val, expr_val)
        #binding_expr_cor.write('{}\t{:3f}\n'.format(binding_id,r_value))
        binding_expr_cor.write('{:2f}\n'.format(r_value,p_value))
#         if binding_id %10000==0:
#             print(binding_id)       
    binding_expr_cor.close()



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--geneid', action = 'store', type = str,dest = 'geneid', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.geneid)
