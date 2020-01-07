import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect

#gained_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/fz_combine_all_features/f1_combined_csv/T-ALL_CTCF_gained_AllFeatures.csv'
#lost_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/fz_combine_all_features/f1_combined_csv/T-ALL_CTCF_lost_AllFeatures.csv'

chrom_size_file = '/nv/vol190/zanglab/zw5j/data/Genome/hg38/hg38_clean.chrom.sizes'
chrom_size_df = pd.read_csv(chrom_size_file,sep='\t',header=None,index_col=0);
chrom_size_df.columns = ['len']

chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']


def main():
    
    #for data in ['Jurkat','A6010']:
    #    for resolution in [5000,20000,100000,1000000]:
    #        for normalization in ['raw','iced']:
    for chrom in chroms:
    
        slurm_dir = 'py1_slurms'
        os.makedirs(slurm_dir,exist_ok=True)
        basename = '{}'.format(chrom)
        outfile = slurm_dir+os.sep+'run_{}.slurm'.format(basename)
        with open(outfile,'w') as slurmout:
            slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
''')
                #######################
                # REMEMBER TO CHANGE THE MEMORY !!!
                #########################
            
            slurmout.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(slurm_dir,basename))
            slurmout.write('time python 1_get_completed_view2M_interaction.py -c {} \n'.format(chrom))
    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-d', '--data', action = 'store', type = str,dest = 'data', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-n', '--normalization', action = 'store', type = str,dest = 'normalization', help = 'input file of', metavar = '<str>')
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
