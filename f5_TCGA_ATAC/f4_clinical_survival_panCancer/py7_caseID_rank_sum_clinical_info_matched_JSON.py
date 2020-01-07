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
import json
    
    
def main_sep_data(clinical_data):
    
    outdir = 'f7_combined_ATAC_sig_diff_score_and_clinical_info_JSON_{}'.format(clinical_data)
    os.makedirs(outdir,exist_ok=True)
    
    # selected cancer-specific ATAC info (overlapped with cancer specific CTCFs)
    atac_sig_file_dir = 'f2_ATAC_sample_irwin_hall'
    cancertypes = ['BRCA','CRC','LUAD','PRAD','PRAD_TissueAdded']
    name_types_cancer_match={'BRCA':'BRCA','CRC':'COAD','LUAD':'LUAD','PRAD':'PRAD','PRAD_TissueAdded':'PRAD'}
    
    for name_type in cancertypes:
        # this is for clinical info for each cancer type
#         clinical_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/15_TCGA_patient_data/clinical_survival/data_TCGA_clinical/clinical.cases_{}.2019-04-08/clinical.tsv'.format(name_types_cancer_match[name_type])
        clinical_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/15_TCGA_patient_data/clinical_survival/data_TCGA_clinical/clinical.project-TCGA-{}.{}.json'.format(name_types_cancer_match[name_type],clinical_data)
        with open(clinical_file) as clinical_inf:
#             clinical_df = pd.read_csv(clinical_inf,sep='\t',index_col=0)
            clinical_list = json.load(clinical_inf)
        clinical_case_id=[i['case_id'] for i in clinical_list]

        # Accessibility and rank-sum (last col) for each caseID
        for binding_type in ['gained','lost']:
#             sig_file = atac_sig_file_dir+os.sep+'{}_{}_irwin_hall_cancertype_ID.csv'.format(name_type,binding_type)   
            sig_file = atac_sig_file_dir+os.sep+'{}_{}_irwin_hall.csv'.format(name_type,binding_type)   
            with open(sig_file) as sig_inf:
                sig_df = pd.read_csv(sig_inf,sep=',',index_col=0)
            sig_df = sig_df.iloc[:,-1:] # keep the rank sum col
            sig_case_id = [i.split('_')[1] for i in sig_df.index] # remove the cancer type prefix
            
            out_df = pd.DataFrame()
            # search by clinical_list rather than sig_df
            # combine clinical info with openness info
            for clinical_case in clinical_list:
                case_id = clinical_case['case_id']
                if case_id in sig_case_id:
                    out_df.loc[case_id,'rank_sum'] = sig_df.loc['{}_{}'.format(name_types_cancer_match[name_type],case_id)]['rank_sum']
                    out_df.loc[case_id,'tumor_stage'] = clinical_case['diagnoses'][0]['tumor_stage']
                    out_df.loc[case_id,'age_at_diagnosis'] = clinical_case['diagnoses'][0]['age_at_diagnosis']
                    out_df.loc[case_id,'vital_status'] = clinical_case['demographic']['vital_status'] 
#                     print(out_df.loc[case_id,'vital_status'])
                    if out_df.loc[case_id,'vital_status']=='Alive':
                        out_df.loc[case_id,'days_to_last_follow_up'] = clinical_case['diagnoses'][0]['days_to_last_follow_up']
                    elif out_df.loc[case_id,'vital_status']=='Dead':
                        out_df.loc[case_id,'days_to_death'] = clinical_case['demographic']['days_to_death']

            
            out_df = out_df.sort_values(by=['rank_sum'],ascending=True)
            out_df.index = ['{}_{}'.format(name_types_cancer_match[name_type],i) for i in out_df.index]
            out_df.to_csv(outdir+os.sep+'{}_{}_irwin_hall_with_clinical_info.csv'.format(name_type,binding_type))

                
def main(infile):

    for clinical_data in ['2019-10-09']:
        main_sep_data(clinical_data)
                


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
