import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.6)
sns.set_style("whitegrid", {'axes.grid' : False})
import CTCF_TALL_modules_new
import scipy
from scipy import stats
import bisect
sns.set_style("ticks")


def write_out_genes(target_genes,deg,flag,flag2,outdir):
    outfile = outdir+os.sep+'{}_{}_allgenes.txt'.format(flag2,flag)
    with open(outfile,'w') as outf:
        outf.write('\n'.join(set(target_genes))+'\n')
    outfile = outdir+os.sep+'{}_{}_upgenes.txt'.format(flag2,flag)
    with open(outfile,'w') as outf:
        outf.write('\n'.join(set(target_genes).intersection(deg))+'\n')

    

def return_intra_domain_genes(ctcf_domain_gene_dic,ctcf_ids):
    intra_domain_genes = set()
    for ctcf_id in ctcf_ids:
        domain_genes = ctcf_domain_gene_dic[str(ctcf_id)]['intra_domain_genes']
        intra_domain_genes = intra_domain_genes.union(domain_genes)
    return intra_domain_genes



def main(outdir):

    outdir='f1_dNOTCH_gainedCTCF_data'
    os.makedirs(outdir,exist_ok=True)
    
    # gained CTCF, intra-domain gene for each CTCF
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    const_df = CTCF_TALL_modules_new.return_constitutive_df()
    ctcf_domain_gene_dic = CTCF_TALL_modules_new.return_union_domain_gene()

    # DEG in different groups
    fcthre=np.log2(1.2)
    pthre=0.05    
    # ==== CUTLL1 up-genes
    cutll1_all,cutll1_upgenes,cutll1_dngenes = CTCF_TALL_modules_new.return_CUTLL1_vs_CD4_deg(fcthre,pthre)  
    tall_all,tall_upgenes,tall_dngenes = CTCF_TALL_modules_new.return_TALL_vs_CD4_deg(fcthre,pthre)
    # ==== NOTCH1 targets
    notch_targets = CTCF_TALL_modules_new.return_Notch1_target_genes()    
    # ==== shCTCF down regulated genes
    shCTCF_all,shCTCF_upgenes,shCTCF_dngenes = CTCF_TALL_modules_new.return_CUTLL1_shCTCF_deg(fcthre,pthre)    
    # ==== GSI down genes
    gsi_all,gsi_upgenes,gsi_dngenes = CTCF_TALL_modules_new.return_CUTLL1_GSI_deg(fcthre,pthre)  
    # ==== GSI washout up genes
    gsi_wo_all,gsi_wo_upgenes,gsi_wo_dngenes = CTCF_TALL_modules_new.return_CUTLL1_GSI_wo_deg(fcthre,pthre)  

    # T-ALL combined feature, get those CTCFs with at least one dynamic NOTCH binding
    combined_df = CTCF_TALL_modules_new.return_cancer_specific_combined_features('T-ALL')
    dynamic_notch_CTCF = combined_df[combined_df['if_intra_domain_dynamic_notch']==1]
    
    # return intra domain genes 
    dynamic_notch_intra_domain_genes = return_intra_domain_genes(ctcf_domain_gene_dic,dynamic_notch_CTCF.index)
    dynamic_notch_gained_CTCF_intra_domain_genes = return_intra_domain_genes(ctcf_domain_gene_dic,dynamic_notch_CTCF.index.intersection(gained_df.index))
    
    
    # overlap with CUTLL1 up-genes
    flag2 = 'cutll1_vs_Tcell'
    target_genes = cutll1_all
    write_out_genes(target_genes,cutll1_upgenes,'all_genes',flag2,outdir)
    
    target_genes = dynamic_notch_intra_domain_genes.intersection(cutll1_all)
    write_out_genes(target_genes,cutll1_upgenes,'dNOTCH_intra_domain',flag2,outdir)
    
    target_genes = dynamic_notch_gained_CTCF_intra_domain_genes.intersection(cutll1_all)
    write_out_genes(target_genes,cutll1_upgenes,'dNOTCH_gainedCTCF_intra_domain',flag2,outdir)
    

    flag2 = 'tall_vs_Tcell'
    target_genes = tall_all
    write_out_genes(target_genes,tall_upgenes,'all_genes',flag2,outdir)
    
    target_genes = dynamic_notch_intra_domain_genes.intersection(tall_all)
    write_out_genes(target_genes,tall_upgenes,'dNOTCH_intra_domain',flag2,outdir)
    
    target_genes = dynamic_notch_gained_CTCF_intra_domain_genes.intersection(tall_all)
    write_out_genes(target_genes,tall_upgenes,'dNOTCH_gainedCTCF_intra_domain',flag2,outdir)
    




    





if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
#     parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.outdir)
