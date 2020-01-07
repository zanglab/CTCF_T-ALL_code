import sys,argparse
import os,glob
import numpy as np
import pandas as pd

import CTCF_TALL_modules_new












def main():


    outdir = 'domain_deg'
    os.makedirs(outdir,exist_ok=True)
    
    deg_info_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f7_panCancer_deg_figs/f3_FeatureCombinaiton_CancerSpecificBinding_local"
    suf_names = ['domain_GT100K_LT1M_log2FC0.585_padj1e-3','domain_GT100K_LT1M_log2FC1_padj1e-5',\
                 'domain_NoLimit_log2FC0.585_padj1e-3','domain_NoLimit_log2FC1_padj1e-5']

    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']
    
    for cancertype in cancertypes:
        for binding_type in ['gained','lost']:
            for suf_name in suf_names:
                info_file = deg_info_dir+os.sep+'f3_cancer_specific_binding_{}/{}_{}.csv'.format(suf_name,cancertype,binding_type)
                info_df = pd.read_csv(info_file,index_col=0)
                
                domain_upgenes = [j for i in info_df['domain_upgenes_names'].dropna().values for j in i.split(',')]
                domain_downgenes = [j for i in info_df['domain_dngenes_names'].dropna().values for j in i.split(',')]
                
                with open(outdir+os.sep+'{}_{}_up_genes_{}.txt'.format(cancertype,binding_type,suf_name),'w') as outf:
                    outf.write('\n'.join(domain_upgenes)+'\n')
                with open(outdir+os.sep+'{}_{}_down_genes_{}.txt'.format(cancertype,binding_type,suf_name),'w') as outf:
                    outf.write('\n'.join(domain_downgenes)+'\n')
                
                






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

