import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import CTCF_TALL_modules_new
import scipy
from scipy import stats


def return_ctcf_intra_domain_deg(gene_set,gained_df,ctcf_domain_gene_dic,outdir,flag):
    
    all_ctcf_genes = set()
    ctcf_deg_df = pd.DataFrame(columns=['intra_domain_{}'.format(flag)])
    for gained_id in gained_df.index:
        domain_genes = ctcf_domain_gene_dic[str(gained_id)]['intra_domain_genes']
        interacted_genes = sorted(gene_set.intersection(domain_genes))
        ctcf_deg_df.loc[gained_id] = ','.join(interacted_genes)
        all_ctcf_genes = all_ctcf_genes.union(interacted_genes)
    
    outfile_base = outdir+os.sep+'gained_CTCF_{}'.format(flag)
    ctcf_deg_df.to_csv(outfile_base+'_df.csv')
    open(outfile_base+'_genes.txt','w').write('\n'.join(all_ctcf_genes)+'\n')
    
    open(outfile_base+'_genes_genome.txt','w').write('\n'.join(gene_set)+'\n')
    




def main():

    outdir = 'f1_intra_domain_deg'
    os.makedirs(outdir,exist_ok=True) 

    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    ctcf_domain_gene_dic = CTCF_TALL_modules_new.return_union_domain_gene()
    fcthre=np.log2(1.2)
    pthre=0.05    
    
    # ==== genome wide genes    
    all_genes = CTCF_TALL_modules_new.return_genome_gene_collections()
    return_ctcf_intra_domain_deg(set(all_genes),gained_df,ctcf_domain_gene_dic,outdir,'all_genes')
    
    # ==== CUTLL1 up-genes
    genes_all,genes_upgenes,genes_dngenes = CTCF_TALL_modules_new.return_CUTLL1_vs_CD4_deg(fcthre,pthre)  
    return_ctcf_intra_domain_deg(set(genes_upgenes),gained_df,ctcf_domain_gene_dic,outdir,'CUTLL1_up')    
    
    # ==== NOTCH1 targets
    notch_targets = CTCF_TALL_modules_new.return_Notch1_target_genes()    
    return_ctcf_intra_domain_deg(set(notch_targets),gained_df,ctcf_domain_gene_dic,outdir,'notch1_target')

    # ==== shCTCF down regulated genes
    shCTCF_all,shCTCF_upgenes,shCTCF_dngenes = CTCF_TALL_modules_new.return_CUTLL1_shCTCF_deg(fcthre,pthre)    
    return_ctcf_intra_domain_deg(set(shCTCF_dngenes),gained_df,ctcf_domain_gene_dic,outdir,'shCTCF_down')
    
    # ==== GSI down genes
    genes_all,genes_upgenes,genes_dngenes = CTCF_TALL_modules_new.return_CUTLL1_GSI_deg(fcthre,pthre)  
    return_ctcf_intra_domain_deg(set(genes_dngenes),gained_df,ctcf_domain_gene_dic,outdir,'GSI_down')
    
    # ==== GSI washout up genes
    genes_all,genes_upgenes,genes_dngenes = CTCF_TALL_modules_new.return_CUTLL1_GSI_wo_deg(fcthre,pthre)  
    return_ctcf_intra_domain_deg(set(genes_upgenes),gained_df,ctcf_domain_gene_dic,outdir,'GSI_wo_up')





 
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

