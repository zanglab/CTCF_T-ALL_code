import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
import scipy
from scipy import stats




def cancer_type_occurrence(peak_occur_in_union_df,cancerType,cancercols,normalcols,outdir):

    df = pd.DataFrame(index = peak_occur_in_union_df.index)
    cancer_df = peak_occur_in_union_df[cancercols];#print(cancer_df);exit(0)
    normal_df = peak_occur_in_union_df[normalcols]  
    controlcols = [i for i in peak_occur_in_union_df.columns if i not in cancercols]
    #control_df = peak_occur_in_union_df[controlcols]
    
    cancer_sum = cancer_df.sum(axis=1)#;print(cancer_sum);exit(0)
    normal_sum = normal_df.sum(axis=1)  
    all_sum = peak_occur_in_union_df.sum(axis=1)
    
    df['cancer_peaks'] = cancer_sum
    df['normal_peaks'] = normal_sum
    df['all_peaks'] = all_sum#;print(df);exit()

    df['cancer_total'] = len(cancercols)
    df['normal_total'] = len(normalcols)
    df['all_total'] = len(peak_occur_in_union_df.columns)#;print(df);exit()
   
    df.to_csv(outdir+os.sep+'{}_binding_occurrence.csv'.format(cancerType))
    


def main():
    
    outdir = 'f1_cancerType_binding_occurrence'
    os.makedirs(outdir,exist_ok=True) 
    
    peak_occur_in_union_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f2_occupancy_on_union_CTCFs/f1_each_data_peak_overlap_union/union_summits_EachDataOverlapInfo.csv'
#     peak_occur_in_union_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f2_occupancy_on_union_CTCFs/f1_each_data_peak_overlap_union/union_summits_EachDataOverlapInfo_head10000.csv'
    with open(peak_occur_in_union_file) as peak_occur_in_union_inf:
        peak_occur_in_union_df = pd.read_csv(peak_occur_in_union_inf,sep="\t",index_col = 0)  

#     cancerType_data = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f0_infile_revised_201808/CancerTypes_datasets.xlsx'
    cancerType_data = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f0_infile/CancerTypes_datasets_201907.xlsx'
    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']

    for cancerType in cancertypes:
        cancertype_df = pd.read_excel(cancerType_data,index_col=0,sheetname=cancerType)
        # cancer and normal GSM IDs
        cancercols = cancertype_df['CancerData'].dropna().values
        normalcols = cancertype_df['NormalData'].dropna().values
        print('\n',cancerType,'\tCancer-len:',len(cancercols),'\tNormal-len:',len(normalcols))
        print('\n',cancercols,'\n\n',normalcols,'\n')
#         exit()
        cancer_type_occurrence(peak_occur_in_union_df,cancerType,cancercols,normalcols,outdir)



       

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of bed fromat, union all the overlapping regions', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
   # parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
