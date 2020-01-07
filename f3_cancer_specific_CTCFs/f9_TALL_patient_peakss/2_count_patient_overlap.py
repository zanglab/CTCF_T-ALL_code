import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
import scipy
from scipy import stats
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['agg.path.chunksize'] = 10000
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=1.1)
#sns.set_style("whitegrid", {'axes.grid' : False})



def main():
    outdir = 'f2_patient_overlap_counting'
    os.makedirs(outdir,exist_ok=True)
    
    all_binding_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3.bed'
    with open(all_binding_file) as all_binding_inf:
        df = pd.read_csv(all_binding_inf,sep='\t',header=None,index_col=3)
    df.columns = ['chr','start','end']
    df.insert(3,'id',df.index)
    #print(df);exit()
    file_pardir = 'patient_overlap_out'    
    for patient in ['PD31','PD9']:
        patient_file = file_pardir+os.sep+'{}_overlap.bed'.format(patient)
        with open(patient_file) as patient_inf:
            patient_df = pd.read_csv(patient_inf,sep='\t',index_col=3,header=None)
        patient_df.columns = ['chr','start','end','{}_peak'.format(patient)]
        df = pd.concat([df,patient_df[['{}_peak'.format(patient)]]],axis=1)
    df['all_patients_peaks'] = df[['PD31_peak','PD9_peak']].sum(axis=1)
    df.to_csv(outdir+os.sep+'union_binding_patient_peak_info.csv',index=False)
                


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
