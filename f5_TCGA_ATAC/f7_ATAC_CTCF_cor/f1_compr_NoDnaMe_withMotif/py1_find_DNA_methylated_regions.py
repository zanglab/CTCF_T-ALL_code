
import os,re,argparse,sys
import bisect
import pandas as pd
from scipy import stats
import association_with_regions
import numpy as np
plus = re.compile('\+')
minus = re.compile('\-')
      


    
def main():    
    
#     outdir = 'f2_union_binding_differential_DNA_methylation'
#     os.makedirs(outdir,exist_ok=True)
    
    methylation_signal_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f8_DNA_methylation/f2_union_binding_methylation_compr/f1_signal_binding_csv'
    jurkat_file = methylation_signal_dir+os.sep+ 'Jurkat_DNAme_unionBinding.csv'
    
    with open(jurkat_file) as cancer_inf:
        cancer_df = pd.read_csv(cancer_inf,index_col=0,sep='\t|,',engine='python',header=None)

    matchness = (cancer_df!=-1)
    cancer_df = cancer_df[matchness]
    
    # DNA methylated regions 
    df = pd.DataFrame(index=cancer_df.index)  
    df['methy_mean_EachSide150bp'] = cancer_df.mean(axis=1)
    df = df[df>15].dropna()
    df.to_csv('Jurkat_DNA_methy_status_GT15.csv')
        
    # DNA methylated regions, with >3 CpGs
    cancer_df_sum = matchness.sum(axis=1)
    cancer_df = cancer_df.loc[cancer_df_sum[cancer_df_sum>3].index]
    df = pd.DataFrame(index=cancer_df.index)  
    df['methy_mean_EachSide150bp'] = cancer_df.mean(axis=1)
    df = df[df>15].dropna()
    df.to_csv('Jurkat_DNA_methy_status_GT15_CpGs_GT3.csv')
 
    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()

