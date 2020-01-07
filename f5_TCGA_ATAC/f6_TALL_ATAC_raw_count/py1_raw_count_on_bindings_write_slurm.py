
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


def write_out_slurm(slurm_dir,bedpe_file,gsm):

    slurmfile = slurm_dir+os.sep+'run_{}.slurm'.format(gsm)
    with open(slurmfile,'w') as slurmout:
        slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
        slurmout.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(slurm_dir,gsm))
        
        py_file = 'get_raw_count_on_regions_readcount.py'
        site_bed_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3.bed'
        outfile_dir = 'f1_TALL_ATAC_raw_count_on_union'
        
        slurmout.write('python {} -s hg38 -f bed -g 0 -i {} -t {} -o {}/{}.txt'.format(py_file,site_bed_file,bedpe_file,outfile_dir,gsm))





def main():
    
    slurm_dir = 'f1_slurm'
    os.makedirs(slurm_dir,exist_ok=True)

    processed_dir='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f13_ATAC_Jurkat_CD4/f1_proprocess'
    
    for gsm in ['GSM1155964','GSM1155965','GSM1155966','GSM1155967','GSM1155968','GSM1155969','GSM2411156','GSM2411157','GSM2411158']:
        bedpe_file = glob.glob('{}/{}/{}_*_PEnoM.bed'.format(processed_dir,gsm,gsm))[0]
        #print(gsm,bedpe_file)
    
        write_out_slurm(slurm_dir,bedpe_file,gsm)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
