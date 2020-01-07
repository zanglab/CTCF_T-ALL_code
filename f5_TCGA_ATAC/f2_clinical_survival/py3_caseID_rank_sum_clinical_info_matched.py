import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

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
from stats_pvalues import irwin_hall_cdf

    
    
def main(infile):
    
    outdir = 'f3_combined_ATAC_sig_diff_score_and_clinical_info'
    os.makedirs(outdir,exist_ok=True)
    
    # selected cancer-specific ATAC info (overlapped with cancer specific CTCFs)
    atac_sig_file_dir = 'f2_ATAC_sample_irwin_hall'
    cancertypes = ['BRCA','CRC','LUAD']#,'PRAD','PRAD_TissueAdded']
    name_types_cancer_match={'BRCA':'BRCA','CRC':'COAD','LUAD':'LUAD','PRAD':'PRAD','PRAD_TissueAdded':'PRAD'}
    
    for name_type in cancertypes:
        # this is for clinical info for each cancer type
        clinical_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/15_TCGA_patient_data/clinical_survival/data_TCGA_clinical/clinical.cases_{}.2019-07-11/clinical.tsv'.format(name_types_cancer_match[name_type])
        clinical_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/15_TCGA_patient_data/clinical_survival/data_TCGA_clinical/clinical.cases_{}.2019-04-08/clinical.tsv'.format(name_types_cancer_match[name_type])
        with open(clinical_file) as clinical_inf:
            clinical_df = pd.read_csv(clinical_inf,sep='\t',index_col=0)
            
        # this is the info for ATAC openness for each caseID
        for binding_type in ['gained','lost']:
            sig_file = atac_sig_file_dir+os.sep+'{}_{}_irwin_hall.csv'.format(name_type,binding_type)   
            with open(sig_file) as sig_inf:
                sig_df = pd.read_csv(sig_inf,sep=',',index_col=0)
            sig_df = sig_df.iloc[:,7:]

            # combine clinical info with openness info
            cancertype_cols = [i for i in sig_df.index if re.search(name_types_cancer_match[name_type],i) and (i.split('_')[1] in clinical_df.index)]
            clinical_tmp = clinical_df.loc[[i.split('_')[1] for i in cancertype_cols]]
            clinical_tmp.index = ['{}_{}'.format(name_types_cancer_match[name_type],i) for i in clinical_tmp.index]
            sig_df_tmp = sig_df.loc[cancertype_cols]
            combined_df = pd.concat([sig_df_tmp,clinical_tmp],axis=1)
            combined_df.to_csv(outdir+os.sep+'{}_{}_irwin_hall_with_clinical_info.csv'.format(name_type,binding_type))

                
             


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
