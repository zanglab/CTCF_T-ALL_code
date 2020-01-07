import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
#sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
#sns.despine(offset=0, trim=True)

import CTCF_TALL_modules_new
            

def find_overlap_by_searching_start_pos(chroms,expand,mutation_df,union_CTCF,outfile):

    outf = open(outfile,'w')
    outf.write('{}\t{}\t{}\t{}\t{}\n'.format('chr','id','mid','strand','mutation'))
    #count_mutation_list = np.zeros(expand*2+1)
    outf2 = open(outfile+'.mutation_event_rate.csv','w')
    for chr in chroms:
        #print(chr)
        mutation_chr = mutation_df[mutation_df['chr']==chr]
        mutation_list = list(mutation_chr['pos'].values)
        mutation_ref2alt = list((mutation_chr['ref']+'>'+mutation_chr['alt']).values) # mutation status/phenotype
            
        ctcf_chr = union_CTCF[union_CTCF['chr']==chr]
        ctcf_list = list(ctcf_chr['mid_position'].values)
        ctcf_index = list(ctcf_chr.index.values)
        ctcf_strand = list(ctcf_chr.motif_strand)
          
        for ctcf_ii in np.arange(len(ctcf_list)):
            # for each CTCF binding site, record the mutation info and mutation rate
            search_phenotype_list = ['']*(expand*2+1)
            mutation_event_rate = np.zeros(expand*2+1)
            ctcf_mid = ctcf_list[ctcf_ii]+1
            ctcf_start = ctcf_mid-expand
            ctcf_end = ctcf_mid+expand
            s = bisect.bisect_left(mutation_list,ctcf_start)
            e = bisect.bisect_right(mutation_list,ctcf_end)
            for jj in np.arange(s,e):
                insert_pos = mutation_list[jj]-ctcf_start
                insert_ele = mutation_ref2alt[jj]
                search_phenotype_list[insert_pos]=insert_ele
                mutation_event_rate[insert_pos]+=1
                #if ctcf_index[ctcf_ii] in gained_df.index:
                #count_mutation_list[insert_pos]+=1
                

            outf.write('{}\t{}\t{}\t{}\t{}\n'.format(chr,ctcf_index[ctcf_ii],ctcf_mid,\
                       ctcf_strand[ctcf_ii],','.join(search_phenotype_list)))
            outf2.write('{}\t{}\n'.format(ctcf_index[ctcf_ii],'\t'.join(map(str,mutation_event_rate))))

    outf.close()
    outf2.close()
    #exit()




def check_rank_list(chroms,check_df,chr_name,col_name):

    for chr in chroms:
        df_chr = check_df[check_df[chr_name]==chr]
        pos_list = df_chr[col_name].values
        for i in np.arange(len(pos_list)-1):
            try:
                assert pos_list[i]<=pos_list[i+1]
            except AssertionError:
                print(chr,pos_list[i-1],pos_list[i],pos_list[i])
                       

def main():

#     indir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/15_TCGA_patient_data/mutation/T_ALL/data/ctcf_sites_output_mutations/'
    indir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f7_ICGC_mutation/T_ALL_new/data/ctcf_sites_output_SNPs'
    outdir = 'f1_mutation_around_CTCF'
    os.makedirs(outdir,exist_ok=True)
    
    chroms = CTCF_TALL_modules_new.return_chroms()
    
    union_CTCF = CTCF_TALL_modules_new.return_occupancy_filtered()
    # check if mid positions are sorted
    check_rank_list(chroms,union_CTCF,'chr','mid_position')
    
    # find mutation around CTCFs  
    suffix='.ctcf.tbl'
    infiles = glob.glob(indir+os.sep+'*{}'.format(suffix))
    for infile in sorted(infiles):
        basename = os.path.basename(infile).split(suffix)[0];print(basename)
        with open(infile) as inf:
            df = pd.read_csv(inf,sep='\t',index_col=None)
        df = df.drop_duplicates()
        check_rank_list(chroms,df,'chr','pos')
        for expand in [9,200]:
            outfile=outdir+os.sep+'{}_mutation_on_CTCF_expand{}.csv'.format(basename,expand);print(expand)
            # get the mutation info on CTCF motif
            find_overlap_by_searching_start_pos(chroms,expand,df,union_CTCF,outfile)
            
            
            
            


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
