import os,sys,argparse
import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
# from GenomeData import *
from operator import itemgetter
#def expand_region(summitlist):
import CTCF_TALL_modules_new 
# import association_with_regions
from collections import Counter


def main():               
    
    outdir='f2_occupancy_score'
    os.makedirs(outdir,exist_ok=True)
    
    collection_df = CTCF_TALL_modules_new.return_collection_df()
    union_df = CTCF_TALL_modules_new.return_union_binding_bed()
    
    overlap_file='f1_each_data_peak_overlap_union/union_summits_EachDataOverlapInfo.csv'
    with open(overlap_file) as inf:
        df = pd.read_csv(inf,sep='\t',index_col=0)

    occupancy_df = pd.concat([union_df,df.sum(axis=1)],axis=1)
    occupancy_df.insert(3,'id',occupancy_df.index)
    occupancy_df.columns = ['chr','start','end','id','score','strand','occupancy_score']
    occupancy_df.to_csv('{}/union_CTCF_occupancy_score.csv'.format(outdir),sep='\t',index=False)
#     print(union_df)

    occupancy_counter = Counter(occupancy_df['occupancy_score'])
    occupancy_count_df = pd.DataFrame()
    for ii in np.arange(1,collection_df.shape[0]+1):
        if ii in occupancy_counter:
            occupancy_count_df.loc[ii,'#bindings']=occupancy_counter[ii]
        else:
            occupancy_count_df.loc[ii,'#bindings']=0
    occupancy_count_df.to_csv('{}/occupancy_count.csv'.format(outdir),sep='\t')


         
if __name__=='__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of bed fromat, union all the overlapping regions', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)

    main()