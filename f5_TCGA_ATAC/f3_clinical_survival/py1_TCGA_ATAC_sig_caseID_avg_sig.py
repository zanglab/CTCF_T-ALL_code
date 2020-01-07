import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
from scipy import stats
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})

import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')



def main(infile):
    
    outdir = 'f1_ATAC_sig_caseID_differential_score'
    os.makedirs(outdir,exist_ok=True)
    
    
    # ATAC identifier
    atac_id_info = pd.read_csv('/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/15_TCGA_patient_data/ATAC_seq/data/TCGA_identifier_mapping.txt',sep='\t',index_col=3)
    caseIDs = sorted(set(atac_id_info.index))
    
    # selected cancer-specific ATAC info (overlapped with cancer specific CTCFs)
    atac_sig_file_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f5_TCGA_ATAC/f1_cancer_specific_CTCF_accessibility/f1_ATAC_sig_abstract/mynorm'
    cancertypes = ['BRCA','CRC','LUAD','PRAD','PRAD_TissueAdded']
#     cancertype_match_names={'BRCA':'BRCA','CRC':'COAD','LUAD':'LUAD','PRAD':'PRAD','PRAD_TissueAdded':'PRAD'}
    
    for cancertype in cancertypes:
        for binding_type in ['gained','lost']:
            print(cancertype)
            sig_file = atac_sig_file_dir+os.sep+'{}_{}_ATAC_sig.bed'.format(cancertype,binding_type)
            #print(sig_file)
            with open(sig_file) as sig_inf:
                sig_df = pd.read_csv(sig_inf,sep='\t')
            # averaged info for each caseID 
            caseID_sig_df = pd.DataFrame(index = sig_df.index) # each row is a binding position
            for caseID in caseIDs:
                caseID_reps = atac_id_info.loc[caseID][['bam_prefix']].values
                if len(caseID_reps.shape)>1:
                    caseID_reps = caseID_reps[:,0]
                caseID_reps = ['_'.join(i.split('-')) for i in caseID_reps]
                # caseID with cancertype+patient-replicate-ID (e.g., ), and matched patient uniq ID (e.g., 00f0f7dd-71de-4e4f-b4b5-df860324f2e8)
                caseID_reps = sig_df.columns.intersection(caseID_reps)#;print(caseID,caseID_reps);exit()
                if len(caseID_reps)>0:
                    cancertype_readed=caseID_reps[0].split('_')[0]
                    caseID_sig_df['{}_{}'.format(cancertype_readed,caseID)] = sig_df[caseID_reps].mean(axis=1)
                
            for insert_col in ['#seqnames','start', 'end', 'name', 'score', 'annotation', 'GC'][::-1]:
                caseID_sig_df.insert(0,insert_col,sig_df[insert_col])
            caseID_sig_df.to_csv(outdir+os.sep+'{}_{}_caseID_ATAC_sig_diff_score.csv'.format(cancertype,binding_type),sep='\t',index=False)
            #exit()
                
             


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
  
    main(args.infile)
