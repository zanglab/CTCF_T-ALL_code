import os,sys,argparse
import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
from GenomeData import *
from operator import itemgetter
#def expand_region(summitlist):
#from reads_count import read_count_on_mapfile


def main():               
    
    ctcf_collection_file='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f0_data_collection_new/GEO_collection_tables/CTCF_peak2000_20190704.xlsx'
    ctcf_collection_df = pd.read_excel(ctcf_collection_file,index_col=0)    
    union_CTCF='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f1_union_all_summits/f1_peak2000_datasets_union_summits/union_summits_fe4_width_150_sorted_natural.bed'
    
    slurm_dir='run_py1_slurm'
    os.makedirs(slurm_dir,exist_ok=True)
    RPKM_dir='f1_RPKM_csv'
    os.makedirs(RPKM_dir,exist_ok=True)
    
    for gsmID in ctcf_collection_df.index:
        bam_file = '/nv/vol190/zanglab/zw5j/projects_data/T_ALL_CTCF_panos/CTCF/bam/{}.bam'.format(gsmID)
        # check existence
        if os.path.isfile('{}/{}.csv'.format(RPKM_dir,gsmID)):
            pass
        else:
             with open('./{}/{}.slurm'.format(slurm_dir,gsmID),'w') as slurmout:
#            with open('./{}/{}_rerun2.slurm'.format(slurm_dir,gsmID),'w') as slurmout:
                slurmout.write('''#!/bin/bash\n
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 6:00:00
#SBATCH -p standard 
#SBATCH -A zanglab
''')
                #######################
                # REMEMBER TO CHANGE THE MEMORY !!!
                #########################
            
                slurmout.write('#SBATCH -o {}/{}.out\n\n'.format(slurm_dir,gsmID))
                slurmout.write('time python count_reads_on_regions.py -i {} -b {} -f bam -s hg38 -o {}/{}.csv'.format(union_CTCF,bam_file,RPKM_dir,gsmID))


         
if __name__=='__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-a', '--infile1', action = 'store', type = str,dest = 'infile1', help = 'input file to be compared/separated', metavar = '<file>')
    #parser.add_argument('-b', '--infile2', action = 'store', type = str,dest = 'infile2', help = 'input file to be compared as basic', metavar = '<file>')
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of bed fromat, union all the overlapping regions', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)

    main()
