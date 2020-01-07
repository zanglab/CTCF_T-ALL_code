
import os,re,argparse,sys
import bisect
import pandas as pd
from scipy import stats
import association_with_regions
import numpy as np
plus = re.compile('\+')
minus = re.compile('\-')
      

def cal_methylation_change_150bp(df,jurkat_df,a6010_df):
    
    ttest_out = stats.ttest_rel(jurkat_df.iloc[:,:],a6010_df.iloc[:,:],axis=1)
    df['stats_EachSide150bp'] = ttest_out[0]
    df['pvalue_EachSide150bp'] = ttest_out[1]
    return df


def cal_methylation_change_motif(df,jurkat_df,a6010_df):
    
    jurkat_motif_df = jurkat_df.iloc[:,151-9-1:151+9]
    a6010_motif_df = a6010_df.iloc[:,151-9-1:151+9]
    
    #ttest_out = stats.ttest_rel(jurkat_motif_df,a6010_motif_df,axis=1)
    #df['stats_motif'] = ttest_out[0]
    #df['pvalue_motif'] = ttest_out[1]
    diff_df = jurkat_motif_df-a6010_motif_df
    df['GT25_motif'] = np.sign(diff_df[diff_df>25]).sum(axis=1)
    df['LT25_motif'] = np.sign(diff_df[diff_df<-25]).sum(axis=1)
    return df


def return_methylation_df_union_binding(methylation_signal_dir,cancer_normal_match_list,cov_thre):

    cancer_file = methylation_signal_dir+os.sep+'{}_cov{}_DNAme_unionBinding.csv'.format(cancer_normal_match_list[0],cov_thre)
    
    normal_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f5_DNA_methylation/T_ALL/f1_signal_binding_csv'
    normal_file = normal_dir+os.sep+'union_binding_A6010_DNAme_UnionBinding.csv'
    
    with open(cancer_file) as cancer_inf, open(normal_file) as normal_inf:
         cancer_df = pd.read_csv(cancer_inf,index_col=0,sep='\t|,',engine='python',header=None)
         normal_df = pd.read_csv(normal_inf,index_col=0,sep='\t|,',engine='python',header=None)
    return cancer_df,normal_df


    
def main():    
    
    methylation_signal_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f5_DNA_methylation/patient/methylation_paired_test/f1_signal_binding_csv'
    
    cancer_normal_match_lists = [['PD31','CD4'],['PD9','CD4'],['PD40','CD4'],['PTBG','CD4']]
    for cancer_normal_match_list in cancer_normal_match_lists:
        for cov_thre in [5,10]:
            cancer_df,normal_df = return_methylation_df_union_binding(methylation_signal_dir,cancer_normal_match_list,cov_thre) 
            # output df
            df = pd.DataFrame(index=cancer_df.index)  
            df = cal_methylation_change_150bp(df,cancer_df,normal_df) 
            df = cal_methylation_change_motif(df,cancer_df,normal_df)
            df.to_csv('f2_paired_test_on_region_level_change_on_motif/union_binding_DNA_methylation_{}_vs_{}_cov{}.csv'.format(cancer_normal_match_list[0],cancer_normal_match_list[1],cov_thre))
        
    
  
    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
