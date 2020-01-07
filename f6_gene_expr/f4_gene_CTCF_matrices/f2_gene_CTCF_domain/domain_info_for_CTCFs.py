import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from GenomeData import *

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
plus = re.compile('\+')
minus = re.compile('\-')
import CTCF_TALL_modules_new


def left_extension(left_point,constitutive_chr):
    # look for a new_left side constitutive CTCF with motif of different direction of left_point 
    new_left_point = left_point
    while constitutive_chr.iloc[new_left_point].motif_strand == constitutive_chr.iloc[left_point].motif_strand:
        if new_left_point==0: # otherwise would out of range for next iteraction
            new_left_point = -1
            break
        new_left_point = new_left_point - 1; 
    return new_left_point
    
def right_extension(right_point,constitutive_chr):
    # look for a new_right side constitutive CTCF with motif of different direction of right_point 
    new_right_point = right_point
    while constitutive_chr.iloc[new_right_point].motif_strand == constitutive_chr.iloc[right_point].motif_strand:
        if new_right_point == len(constitutive_chr.index)-1:
            new_right_point = len(constitutive_chr.index)
            break
        new_right_point = new_right_point + 1   
    return new_right_point

def nearby_consitiutive_binding_both_dir(position,constitutive_chr, constitutive_start, constitutive_end):
    # return nearest constitutive CTCF binding
    left_point = bisect.bisect_left(constitutive_start,position) -1
    right_point = bisect.bisect_right(constitutive_end,position)
    # initiate the left/right bindings
    if left_point>=0 and right_point<len(constitutive_chr.index):
        left_binding = constitutive_chr.iloc[left_point]
        right_binding = constitutive_chr.iloc[right_point]
        # if is in different directions, return
        if left_binding.motif_strand != right_binding.motif_strand:
            return left_binding,right_binding,None  
        # search matched binding respectively    
        else:
            # initiate domain len
            newleft_domain_len,newright_domain_len = 999999999,999999999
            strand = left_binding.motif_strand           
            new_left_point = left_extension(left_point,constitutive_chr)
            new_right_point = right_extension(right_point,constitutive_chr)             
            # if have valid new left/right point
            if new_left_point>=0:        
                left_binding1,right_binding1 = constitutive_chr.iloc[new_left_point],constitutive_chr.iloc[right_point]
                newleft_domain_len = int(right_binding1.motif_end) - int(left_binding1.motif_start)
            if new_right_point<len(constitutive_chr.index):
                left_binding2,right_binding2 = constitutive_chr.iloc[left_point],constitutive_chr.iloc[new_right_point]
                newright_domain_len = int(right_binding2.motif_end) - int(left_binding2.motif_start)
            # compr and return smaller domain
            if newleft_domain_len < newright_domain_len:# and newleft_domain_len<5000000:
                return left_binding1,right_binding1,None
            elif newleft_domain_len > newright_domain_len:# and newright_domain_len<5000000:
                return left_binding2,right_binding2,None
                # if None of both expansion work, return -1
            else:
                return None,None,None
    else:
       return None,None,'ChrEnd'       




def main():
    
    
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', \
    'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    # ==== constitutive bindings with motif
    constitutive_motif_binding_file = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_constitutive_0.8.csv"
    constitutive_motif_binding = pd.read_csv(constitutive_motif_binding_file,sep = ',',index_col=3)
    constitutive_motif_binding = constitutive_motif_binding[constitutive_motif_binding['motif_start']!="N"]

    # ==== output file
    outfile = open('all_CTCF_domainInfo.csv','w')
    outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('chr','middle','id','domain_left','domain_right','domain_len'))
    
    # ==== union bindings with occupancy score ≥ 3
    union_df = CTCF_TALL_modules_new.return_occupancy_filtered()

#     chroms = ['chr22']
    for chr in chroms:
        union_df_chr = union_df[union_df['chr']==chr]
        constitutive_chr = constitutive_motif_binding[constitutive_motif_binding['chr']==chr]
        constitutive_start = constitutive_chr['motif_start'].astype(int).tolist()
        constitutive_end = constitutive_chr['motif_end'].astype(int).tolist()

        for ctcf_id in union_df_chr.index:
            ctcf_binding = union_df_chr.loc[ctcf_id]
            # binding_chr = ctcf_binding.chr
            middle = ctcf_binding.mid_position
        
            left_binding,right_binding,chrEndInfo = nearby_consitiutive_binding_both_dir(middle,constitutive_chr,constitutive_start,constitutive_end)
            domain_left,domain_right,domain_len = 'NA','NA','NA'   
            if left_binding is not None:
                domain_left = int(left_binding.motif_start)
                domain_right = int(right_binding.motif_end)
                domain_len = int(right_binding.motif_end) - int(left_binding.motif_start)
            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr,middle,ctcf_id,domain_left,domain_right,domain_len))   
            if ctcf_id%1000==0:
                print(chr,ctcf_id)
    outfile.close()

    
            

 


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--geneid', action = 'store', type = str,dest = 'geneid', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
