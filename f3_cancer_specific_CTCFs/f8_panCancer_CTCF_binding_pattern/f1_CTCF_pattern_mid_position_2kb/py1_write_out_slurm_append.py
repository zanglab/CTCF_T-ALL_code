'''
Created on 04/15/2018
@authors: Chongzhi Zang, Zhenjia Wang<zw5j@virginia.edu>

This file is used to generate pattern collected from input tagsfile of bed/bam format 
for each start position width-expanded-site in input sitefile
'''

import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
#import association_with_regions
from get_reads_positions import reads_positions
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"

import re,bisect
plus = re.compile('\+')
minus = re.compile('\-')


def write_out_slurm(slurm_dir,bamfile,cancerType,flag,gsm_name):

    slurmfile = slurm_dir+os.sep+'run_{}_{}_{}.slurm'.format(cancerType,flag,gsm_name)
    with open(slurmfile,'w') as slurmout:
        slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=60000
#SBATCH -t 12:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
        slurmout.write('#SBATCH -o {}/slurm_{}_{}_{}.out\n\n'.format(slurm_dir,cancerType,flag,gsm_name))
        py_file = 'get_pattern_near_site_readcount_write_out.py'
        site_bed_file = 'union_binding_occupancy_score_GT3_mid_position.bed'
        outfile_dir = 'f1_ctcf_binding_pattern_csv'
        
        slurmout.write('python {} -s hg38 -f bam -m -w 1000 -b 200 -i {} -t {} -o {}/{}_{}_{}.txt'.format(py_file,site_bed_file,bamfile,outfile_dir,cancerType,flag,gsm_name))



   
def main():

    slurm_dir = 'f0_slurm_append'
    os.makedirs(slurm_dir,exist_ok=True)

    bam_file_dir = '/nv/vol190/zanglab/zw5j/projects_data/T_ALL_CTCF_panos/CTCF/bam'
    
    cancerType_data = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f0_infile/CancerTypes_datasets_201907.xlsx'
    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']

    for cancerType in cancertypes:
        cancertype_df = pd.read_excel(cancerType_data,index_col=0,sheetname=cancerType)
        # cancer and normal GSM IDs
        cancercols = cancertype_df['CancerData'].dropna().values
        normalcols = cancertype_df['NormalData'].dropna().values
        print('\n## ====',cancerType,'Cancer-len:',len(cancercols),'\tNormal-len:',len(normalcols))
#         print('\n',cancercols,'\n',normalcols,'\n')

        for bamID in cancercols:
            bam_file = bam_file_dir+os.sep+'{}.bam'.format(bamID)
#             print(bam_file,os.path.isfile(bam_file))
#             print('python get_pattern_near_site_readcount_write_out.py -s hg38 -f bam -m -w 1000 -b 200 -i union_binding_occupancy_score_GT3_mid_position.bed -t {} -o f1_ctcf_binding_pattern_csv/{}_{}_{}.txt'.format(bam_file,cancerType,'CancerCol',bamID))
            if not os.path.isfile('f1_ctcf_binding_pattern_csv/{}_cancer_{}.txt'.format(cancerType,bamID)):
                write_out_slurm(slurm_dir,bam_file,cancerType,'cancer',bamID)
                
        for bamID in normalcols:
            bam_file = bam_file_dir+os.sep+'{}.bam'.format(bamID)
#             print(bam_file,os.path.isfile(bam_file))
#             print('python get_pattern_near_site_readcount_write_out.py -s hg38 -f bam -m -w 1000 -b 200 -i union_binding_occupancy_score_GT3_mid_position.bed -t {} -o f1_ctcf_binding_pattern_csv/{}_{}_{}.txt'.format(bam_file,cancerType,'NormCol',bamID))
            if not os.path.isfile('f1_ctcf_binding_pattern_csv/{}_normal_{}.txt'.format(cancerType,bamID)):
                write_out_slurm(slurm_dir,bam_file,cancerType,'normal',bamID)







if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
