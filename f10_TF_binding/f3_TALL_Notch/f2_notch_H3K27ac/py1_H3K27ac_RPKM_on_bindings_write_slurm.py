
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


def write_out_slurm(slurm_dir,bamfile,celltype,hm,prename):

    slurmfile = slurm_dir+os.sep+'run_{}_{}_{}.slurm'.format(celltype,hm,prename)
    with open(slurmfile,'w') as slurmout:
        slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
        slurmout.write('#SBATCH -o {}/slurm_{}_{}_{}.out\n\n'.format(slurm_dir,celltype,hm,prename))
        py_file = 'get_RPKM_on_regions_readcount.py'
        site_bed_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3_mid_position.bed'
        outfile_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f10_TF_binding/f3_TALL_Notch/f2_notch_H3K27ac/f1_HM_RPKM'
        
        slurmout.write('python {} -s hg38 -f bam -m -e 10000 -i {} -t {} -o {}/{}_{}_{}.txt'.format(py_file,site_bed_file,bamfile,outfile_dir,celltype,hm,prename))



   

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
    
    slurm_dir = 'f1_slurm'
    os.makedirs(slurm_dir,exist_ok=True)

    celltypes=['CD4','JURKAT']
    for celltype in celltypes:
        df = pd.read_excel('/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f11_histone_modification/data_files/CTCF_pancancer_histone_modification.xlsx',sheet_name=celltype)
        for hm in ['H3K27ac']:
            prenames = df[hm].dropna()
            for prename in prenames:
                bamfile,prename = return_hm_bam_file(prename)
                print(celltype,os.path.basename(bamfile))
                #if not os.path.isfile('{}/{}_{}_{}_{}.txt'.format('HM_pattern',cancertype,celltype,hm,prename)):
                write_out_slurm(slurm_dir,bamfile,celltype,hm,prename)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
