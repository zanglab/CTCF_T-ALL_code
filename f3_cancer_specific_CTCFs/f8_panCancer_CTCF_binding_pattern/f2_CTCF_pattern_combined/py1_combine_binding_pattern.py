import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
#import association_with_regions
from get_reads_positions import reads_positions
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"



   
def main():

    outdir = 'f1_combined_cancerType_CTCF_pattern_csv'
    os.makedirs(outdir,exist_ok=True)
    
    # averaged dataframe for each cancertype and each datatype
    csv_file_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f8_panCancer_CTCF_binding_pattern/f1_CTCF_pattern_mid_position_2kb/f1_ctcf_binding_pattern_csv'
    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']
    for cancertype in cancertypes:
        for datatype in ['cancer','normal']:
            csv_files = glob.glob('{}/{}_{}_*'.format(csv_file_dir,cancertype,datatype)) 
            df = pd.DataFrame()
            for csv_file in csv_files:
                print(cancertype,datatype,len(csv_files),os.path.basename(csv_file))#;exit() 
                with open(csv_file) as csv_inf:
                    csv_df = pd.read_csv(csv_inf,sep="\t",index_col=0)
                if df.shape[0]==0:
                    df = csv_df
                else:
                    df = df+csv_df
            df = df/len(csv_files) 
            df = df.round(2)     
            df.to_csv(outdir+os.sep+'{}_{}_avg_combined.csv'.format(cancertype,datatype))

    # averaged dataframe for T-ALL Jurkat and CUTLL1
    csv_file_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f8_panCancer_CTCF_binding_pattern/f1_CTCF_pattern_mid_position_2kb/f1_ctcf_binding_pattern_csv'
    if 1:
        for prename in ['T-ALL_cancer_JURKAT','T-ALL_cancer_CUTLL1']:
            csv_files = glob.glob('{}/{}_*'.format(csv_file_dir,prename)) 
            df = pd.DataFrame()
            for csv_file in csv_files:
                print(prename,len(csv_files),os.path.basename(csv_file))#;exit() 
                with open(csv_file) as csv_inf:
                    csv_df = pd.read_csv(csv_inf,sep="\t",index_col=0)
                if df.shape[0]==0:
                    df = csv_df
                else:
                    df = df+csv_df
            df = df/len(csv_files) 
            df = df.round(2)     
            df.to_csv(outdir+os.sep+'{}_avg_combined.csv'.format(prename))


    # averaged dataframe for T-ALL patient samples
    csv_file_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f8_panCancer_CTCF_binding_pattern/f1_CTCF_pattern_mid_position_2kb/f2_CTCF_pattern_csv_patient'
    if 1:
        for prename in ['PD9','PD31']:
            csv_files = glob.glob('{}/{}_*'.format(csv_file_dir,prename)) 
            df = pd.DataFrame()
            for csv_file in csv_files:
                print(prename,len(csv_files),os.path.basename(csv_file))#;exit() 
                with open(csv_file) as csv_inf:
                    csv_df = pd.read_csv(csv_inf,sep="\t",index_col=0)
                if df.shape[0]==0:
                    df = csv_df
                else:
                    df = df+csv_df
            df = df/len(csv_files) 
            df = df.round(2)     
            df.to_csv(outdir+os.sep+'{}_avg_combined.csv'.format(prename))




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
