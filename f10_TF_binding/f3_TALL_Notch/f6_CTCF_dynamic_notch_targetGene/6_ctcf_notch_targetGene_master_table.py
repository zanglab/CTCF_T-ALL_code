import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import association_with_genes
import association_with_regions
import re,bisect
import CTCF_TALL_modules_new



def return_tf_gene_pairwise_dis(master_df,col_name,mark):

    tf_gbstart = master_df[col_name]-master_df['gbStart']
    tf_gbend = master_df[col_name]-master_df['gbEnd']
    if_gb_overlap_sign = tf_gbstart*tf_gbend
    gb_Nonoverlapped = if_gb_overlap_sign[if_gb_overlap_sign>0].index
    master_df.loc[gb_Nonoverlapped,'dis_{}_gene'.format(mark)] = np.min((np.abs(tf_gbstart[gb_Nonoverlapped]),np.abs(tf_gbend[gb_Nonoverlapped])),axis=0)
    master_df['dis_{}_promoter'.format(mark)] = np.abs(master_df[col_name]-master_df['promoter'])
    
    return master_df


def main():

    outdir = 'f6_ctcf_notch_gene_master'
    os.makedirs(outdir,exist_ok=True) 

    fcthre=np.log2(1.2)
    pthre=0.05    
    # ==== CUTLL1 up-genes
    cutll1_all,cutll1_upgenes,cutll1_dngenes = CTCF_TALL_modules_new.return_CUTLL1_vs_CD4_deg(fcthre,pthre)  
    # ==== NOTCH1 targets
    notch_targets = CTCF_TALL_modules_new.return_Notch1_target_genes()    
    # ==== shCTCF down regulated genes
    shCTCF_all,shCTCF_upgenes,shCTCF_dngenes = CTCF_TALL_modules_new.return_CUTLL1_shCTCF_deg(fcthre,pthre)    
    # ==== GSI down genes
    gsi_all,gsi_upgenes,gsi_dngenes = CTCF_TALL_modules_new.return_CUTLL1_GSI_deg(fcthre,pthre)  
    # ==== GSI washout up genes
    gsi_wo_all,gsi_wo_upgenes,gsi_wo_dngenes = CTCF_TALL_modules_new.return_CUTLL1_GSI_wo_deg(fcthre,pthre)  
    
    # gene ucsc format
    ucsc_df = association_with_genes.return_ucsc_df('hg38')
    # intra-domain genes for each CTCF
    ctcf_domain_gene_dic = CTCF_TALL_modules_new.return_union_domain_gene()
    combined_df = CTCF_TALL_modules_new.return_cancer_specific_combined_features('T-ALL')
    
    # ctcf with intra-notch1 binding
    z_score_file = 'f4_interaction_compr_zscores/raw_T-ALL_gained_zscore_delta_change_filtered.txt'
    z_score_df = pd.read_csv(z_score_file,sep='\t',index_col=0)

    # for each CTCF, get intra-domain gene and gene expression pattern
    master_df = pd.DataFrame()
    master_df_ii=0
    for gained_id in z_score_df.index:
        domain_genes = ctcf_domain_gene_dic[str(gained_id)]['intra_domain_genes']
        cutll1_deg = sorted(set(cutll1_upgenes).intersection(domain_genes))
        for gene in cutll1_deg:
            master_df.loc[master_df_ii,'chr'] = combined_df.loc[gained_id,'chr']
            master_df.loc[master_df_ii,'ctcf_id'] = gained_id
            master_df.loc[master_df_ii,'ctcf_pos'] = combined_df.loc[gained_id,'mid_position']
            master_df.loc[master_df_ii,'ctcf_annotation'] = combined_df.loc[gained_id,'annotation']
            master_df.loc[master_df_ii,'hicor_notch_pos'] = z_score_df.loc[gained_id,'notch_pos']
            master_df.loc[master_df_ii,'gene_id'] = gene
            master_df.loc[master_df_ii,'gene_strand']  = ucsc_df.loc[gene].strand
            master_df.loc[master_df_ii,'gbStart']  = ucsc_df.loc[gene].gbStart
            master_df.loc[master_df_ii,'gbEnd']  = ucsc_df.loc[gene].gbEnd
            # add promoter
            if ucsc_df.loc[gene].strand=='+':
                master_df.loc[master_df_ii,'promoter']  = ucsc_df.loc[gene].gbStart
            else:
                master_df.loc[master_df_ii,'promoter']  = ucsc_df.loc[gene].gbEnd
            # add DEG info
            if gene in notch_targets:
                master_df.loc[master_df_ii,'notch_target']  = 1
            if gene in shCTCF_dngenes:
                master_df.loc[master_df_ii,'shCTCF_dngenes']  = 1
            if gene in gsi_dngenes:
                master_df.loc[master_df_ii,'gsi_dngenes']  = 1
            if gene in gsi_wo_upgenes:
                master_df.loc[master_df_ii,'gsi_wo_upgenes']  = 1
            master_df_ii+=1
    
    master_df = return_tf_gene_pairwise_dis(master_df,'ctcf_pos','ctcf')
    master_df = return_tf_gene_pairwise_dis(master_df,'hicor_notch_pos','notch')
    master_df['dis_ctcf_notch'] = np.abs(master_df['ctcf_pos']-master_df['hicor_notch_pos'])
    master_df['gene_len'] = np.abs(master_df['gbStart']-master_df['gbEnd'])
    
    master_df = master_df.fillna(0)
    master_df.to_csv(outdir+os.sep+'T_ALL_gained_CTCF_NOTCH_gene_master.csv')
            
            


 
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

