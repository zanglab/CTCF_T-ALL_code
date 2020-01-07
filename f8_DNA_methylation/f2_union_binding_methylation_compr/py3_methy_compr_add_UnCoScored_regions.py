
import os,re,argparse,sys
import bisect
import pandas as pd
from scipy import stats
import association_with_regions
import numpy as np
plus = re.compile('\+')
minus = re.compile('\-')
      

def cal_methylation_change_150bp(df,cancer_df,normal_df):
    
    df['differential_methy_EachSide150bp'] = cancer_df.mean(axis=1)-normal_df.mean(axis=1)
    df['DE_count_EachSide150bp'] = cancer_df.notnull().sum(axis=1)
    return df


def cal_methylation_change_motif(df,cancer_df,normal_df):
    
    cancer_df = cancer_df.iloc[:,151-9-1:151+9]
    normal_df = normal_df.iloc[:,151-9-1:151+9]
    
    df['differential_methy_motif'] = cancer_df.mean(axis=1)-normal_df.mean(axis=1)
    df['DE_count_motif'] = cancer_df.notnull().sum(axis=1)
    return df



def return_methylation_df_union_binding(methylation_signal_dir,cancer_normal_match_list):
    cancer_file = methylation_signal_dir+os.sep+'{}_DNAme_unionBinding.csv'.format(cancer_normal_match_list[0])
    normal_file = methylation_signal_dir+os.sep+'{}_DNAme_unionBinding.csv'.format(cancer_normal_match_list[1])
#     print(os.path.isfile(cancer_file),os.path.isfile(normal_file))
    with open(cancer_file) as cancer_inf, open(normal_file) as normal_inf:
         cancer_df = pd.read_csv(cancer_inf,index_col=0,sep='\t|,',engine='python',header=None)
         normal_df = pd.read_csv(normal_inf,index_col=0,sep='\t|,',engine='python',header=None)
    return cancer_df,normal_df


    
def main():    
    
    outdir = 'f3_union_binding_differential_DNA_methylation_append'
    os.makedirs(outdir,exist_ok=True)
    
    methylation_signal_dir = 'f1_signal_binding_csv'
    cancer_normal_match_lists=\
    [['A549_cov5','lung_cov5_WGBS'],['A549_cov5','lung_cov5'],\
    ['HCT116_cov5','sigmoid_colon_cov5'],['HCT116_GSM1465024_cov5_WGBS','sigmoid_colon_cov5'],\
    ['HCT116_GSM257964x_cov5_WGBS','sigmoid_colon_cov5'],\
    ['MCF7_cov5','breast_cov5'],['CUTLL1','A6010'],['Jurkat','A6010'],\
    ['PD31_cov5','A6010'],['PD9_cov5','A6010'],\
    ['LNCaP_cov5','prostate_gland_cov5']]
    

#     cancer_normal_match_lists = [['LNCaP_cov5','prostate_gland_cov5']]



    for cancer_normal_match_list in cancer_normal_match_lists:
        cancer_df,normal_df = return_methylation_df_union_binding(methylation_signal_dir,cancer_normal_match_list) 
        matchness = (cancer_df!=-1)&(normal_df!=-1)
        matchness_all = (cancer_df!=-1)|(normal_df!=-1)
        cancer_df = cancer_df[matchness]
        normal_df = normal_df[matchness]
        
        df = pd.DataFrame(index=cancer_df.index)  
        df['UnionScored_count_EachSide150bp'] = matchness_all.sum(axis=1)
        df = cal_methylation_change_150bp(df,cancer_df,normal_df) 
        df = cal_methylation_change_motif(df,cancer_df,normal_df)
        df.to_csv('{}/{}_vs_{}.csv'.format(outdir,cancer_normal_match_list[0],cancer_normal_match_list[1]))
        
    
 
    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()

