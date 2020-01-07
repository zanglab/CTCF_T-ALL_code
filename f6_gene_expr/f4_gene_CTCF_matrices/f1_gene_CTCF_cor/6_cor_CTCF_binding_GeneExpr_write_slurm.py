import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
import association_with_genes
import association_with_regions
import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
import read_write_file

def write_out_slurm(gene_ids):
    
    with open('f6_slurm_out/{}.slurm'.format(gene_ids[0]),'w') as slurm_out:
        slurm_out.write('''#!/bin/bash\n
#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
''')
                #######################
                # REMEMBER TO CHANGE THE MEMORY !!!
                #########################
            
        slurm_out.write('#SBATCH -o {}.out\n\n'.format(gene_ids[0]))
        for gene_id in gene_ids:
            slurm_out.write('time python ../7_cor_CTCF_binding_GeneExpr.py -i {} \n'.format(gene_id))




def main():
    outdir='f6_slurm_out'
    os.makedirs(outdir,exist_ok=True)
    
    ucsc_df = association_with_genes.return_ucsc_df('hg38')
    gene_expr_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f2_gene_expr_TPM/f5_gene_expr_combined/geneExpr_replicates_combined_TPM_sqrt.csv"
    gene_expr = pd.read_csv(gene_expr_file,index_col = 0)    
    
    gene_collections = ucsc_df.index.intersection(gene_expr.index)
    step=6
    
    for gene_id_pos in np.arange(0,len(gene_collections),step):
        gene_ids = gene_collections[gene_id_pos:gene_id_pos+step]
        csv_file = 'f7_cor_CTCF_binding_geneExpr/CTCF_gene_cor_csv/{}_cor.csv'.format(gene_ids[0])
        #if os.path.isfile(csv_file):
        #    csv_lines = read_write_file.get_lines(csv_file)
        #    if csv_lines != 564030:
        #        print(gene_id,'not full')#;exit(0)
        #        write_out_slurm(gene_id)
        if not os.path.isfile(csv_file):
            #print(gene_id,'not file')
            write_out_slurm(gene_ids)

        #print(gene_id)



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
