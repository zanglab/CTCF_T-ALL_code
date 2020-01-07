import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
sns.despine(offset=0, trim=True)
# import CTCF_TALL_modules
import twobitreader
import Bio
from Bio.Seq import Seq
           


def get_mutation_info(infile):
    
    with open(infile) as inf:
        df = pd.read_csv(inf,sep='\t',index_col=1)
    info_len = len(df.loc[df.index[0],'mutation'].split(','))  
    info_df = pd.DataFrame(index = df.index,columns = np.arange(info_len)) 
    for ii in df.index:
        info = df.loc[ii,'mutation'].split(',')   
        ## reverse info on the minus strand
        if df.loc[ii,'strand']=='-':
           info = info[::-1]
        info_df.loc[ii]=info
    avg_mutation_percentage = [len(i)/(3*info_df.shape[0]) for i in info_df.sum()]    
    return info_df,avg_mutation_percentage
        


def if_all_ATCG(sequence):
    countA = sequence.count('A')
    countT = sequence.count('T')
    countC = sequence.count('C')
    countG = sequence.count('G')
    if countA+countT+countC+countG == len(sequence):
        return True
    else:
        return False
    

def get_altered_sequencing(sequence,mutation_info,mutation_info_len):
    # get the ref->alt sequencing
    alt_sequencing = sequence
    for mutation_pos in np.arange(len(mutation_info_len)): # scan the mutation len
        if mutation_info_len[mutation_pos]!=0:
            assert alt_sequencing[mutation_pos] == mutation_info[mutation_pos][0]
            alt_sequencing = alt_sequencing[:mutation_pos]+mutation_info[mutation_pos][-1]+alt_sequencing[mutation_pos+1:]
    return alt_sequencing

def return_matching_score(log2LikelihoodRatio,seq):
#     print(seq)
    score = sum([log2LikelihoodRatio.loc[i,seq[i][0].upper()] for i in np.arange(len(seq)) if seq[i][0] in ['A','C','G','T']])
    rev_seq = seq.reverse_complement()
    rev_score = sum([log2LikelihoodRatio.loc[i,rev_seq[i][0].upper()] for i in np.arange(len(rev_seq))if rev_seq[i][0] in ['A','C','G','T']])
    return score, rev_score


def main():

    outdir = 'f6_mutation_matrix_score'
    os.makedirs(outdir,exist_ok=True)
    indir = 'f2_mutation_on_cancer_specific_CTCF'
    
    # CTCF position weight matrix
    pwm_file ='/nv/vol190/zanglab/shared/Motif/pwm/jaspar_vertebrates/CTCF_MA0139.1.txt' # ACGT
    pwm = pd.read_csv(pwm_file,header=None,sep='\t')
    pwm.columns=['A','C','G','T']
    
    # this is only tested bg_var from /nv/vol190/zanglab/shared/Motif/sites/hg38_fimo_jarspar/raw_results/CTCF/fimo.txt
    bg_var=0.025
    bg=[0.25+bg_var,0.25-bg_var,0.25-bg_var,0.25+bg_var]
    log2LikelihoodRatio = np.log2((pwm+0.0000001)/bg)
    # hg38 sequencing
    genome = twobitreader.TwoBitFile('/nv/vol190/zanglab/zw5j/work2017/fusion_analysis/111_fusion_append/4_CTCFpairing_GCskew_RLoop_RNAprotein_GTExFusion/f0_infiles/hg38.2bit')

    cancertype_mutation_matchness = {'BRCA':'BRCA','CRC':'COAD','LUAD':'LUAD','PRAD':'PRAD','AML':'AML','PRAD_TissueAdded':'PRAD'}
    cancertypes=['BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']
    
    for celltype in cancertypes:
        for binding_type in ['gained','lost']:
            mutation_info_file ='{}/{}_mutation_on_{}_{}_expand9.csv'.format(indir,cancertype_mutation_matchness[celltype],celltype,binding_type)
            df = pd.read_csv(mutation_info_file,sep='\t',index_col=1,low_memory=False)
            # for each binding position with motif, get the DNA sequence
            matching_score_df = pd.DataFrame()
            for id in df.index:
                chr = df.loc[id,'chr']
                mid = df.loc[id,'mid']
                strand = df.loc[id,'strand']
                mutation_info = df.loc[id,'mutation'].split(',')
                mutation_info_len = [len(i) for i in mutation_info]
                # check if there is a mutation
                matching_score_df.loc[id,'chr'] = chr
                matching_score_df.loc[id,'mid'] = mid
                matching_score_df.loc[id,'strand'] = strand
                
                if sum(mutation_info_len) !=0:
                    sequence = genome[chr][mid-9-1:mid+9];sequence
                    sequence = Seq(sequence).upper()
                    assert if_all_ATCG(sequence)
                    alt_sequence = get_altered_sequencing(sequence,mutation_info,mutation_info_len)
                    seq_score,rev_seq_score = return_matching_score(log2LikelihoodRatio,sequence)
                    alt_seq_score,rev_alt_seq_score = return_matching_score(log2LikelihoodRatio,alt_sequence)
                    #matching_score_df[id,'sequence'] = sequence
                   # matching_score_df[id,'alt_sequence'] = alt_sequence
                    matching_score_df.loc[id,'seq_score'] = seq_score
                    matching_score_df.loc[id,'alt_seq_score'] = alt_seq_score
                    matching_score_df.loc[id,'rev_seq_score'] = rev_seq_score
                    matching_score_df.loc[id,'rev_alt_seq_score'] = rev_alt_seq_score
                else:
                    matching_score_df.loc[id,'seq_score'] = 0
                    matching_score_df.loc[id,'alt_seq_score'] = 0
                    matching_score_df.loc[id,'rev_seq_score'] = 0
                    matching_score_df.loc[id,'rev_alt_seq_score'] = 0
                   
            matching_score_df.to_csv('{}/{}_{}.csv'.format(outdir,celltype,binding_type))#;exit()
                    
                        

            


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
