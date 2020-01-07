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

cancer_type_matched = {'T-ALL':['CD4','JURKAT','CUTLL1','PD9','PD31'],'BRCA':['Breast_normal','BRCA'],'LUAD':['Lung_normal','LUAD',],'CRC':['Colon_normal','CRC'],'PRAD':['CRC'],'AML':['AML']}


def write_out_slurm(slurm_dir,bamfile,cancertype,celltype,hm,prename):

    slurmfile = slurm_dir+os.sep+'run_{}_{}_{}_{}.slurm'.format(cancertype,celltype,hm,prename)
    with open(slurmfile,'w') as slurmout:
        slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
        slurmout.write('#SBATCH -o {}/slurm_{}_{}_{}_{}.out\n\n'.format(slurm_dir,cancertype,celltype,hm,prename))
        py_file = 'get_pattern_near_site_readcount_write_out.py'
        site_bed_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3_mid_position.bed'
        outfile_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f9_histone_modification/f1_binding_pattern_on_union/f1_HM_pattern_csv'
        
        slurmout.write('time python {} -s hg38 -f bam -m -w 1000 -b 200 -i {} -t {} -o {}/{}_{}_{}_{}.txt\n\n'.format(py_file,site_bed_file,bamfile,outfile_dir,cancertype,celltype,hm,prename))



   

def return_hm_bam_file(prename):

    # get the bam file for each prename   
    try:
        prename = int(prename)
    except:
        pass
    file1 = '/nv/vol190/zanglab/cz3d/cistrome_data/bam/{}_treat_rep1.bam'.format(prename)#;print(file1)
    file2 = '/nv/vol190/zanglab/zw5j/projects_data/T_ALL_CTCF_panos/Public/histone_modification/{}.bam'.format(prename)
    file3 = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/7_panCancer_DNAmethylation/f2_callpeak/f2_histone_callpeak_chilin/{}/attic/{}_treatment.bam'.format(prename,prename)
    for isfile in [file1,file2,file3]:
        if os.path.isfile(isfile):
            #print(os.path.basename(isfile))
            return isfile,prename
    return False,prename


def main():
    
    slurm_dir = 'f0_slurm'
    csv_dir = 'f1_HM_pattern_csv'
    os.makedirs(slurm_dir,exist_ok=True)
    os.makedirs(csv_dir,exist_ok=True)

    for cancertype in cancer_type_matched.keys():
        for celltype in cancer_type_matched[cancertype]:
            df = pd.read_excel('/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f9_histone_modification/data_file/CTCF_pancancer_histone_modification_201907.xlsx',sheet_name=celltype)
            #print(celltype,df)#;exit()
            for hm in ['H3K27ac','H3K27me3','H3K4me1']:
                prenames = df[hm].dropna()
                for prename in prenames:
                    bamfile,prename = return_hm_bam_file(prename)
                    #print(cancertype,celltype,hm,os.path.basename(bamfile))
                    if not os.path.isfile('{}/{}_{}_{}_{}.txt'.format(csv_dir,cancertype,celltype,hm,prename)):
                        write_out_slurm(slurm_dir,bamfile,cancertype,celltype,hm,prename)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
