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
import association_with_genes
import association_with_regions
import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
#sys.path.insert(0,os.path.abspath('modules'))
#import cancer_specific_binding_identification
import CTCF_TALL_modules_new


def abstract_binding_pattern_df(df,cancertype,prenames,outdir,flag):

    combined_file_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f9_histone_modification/f2_combine_samples/f1_combined_cancerType_csv'
    for prename in prenames:
        for hm in ['H3K27ac','H3K27me3','H3K4me1']:
            infile = combined_file_dir+os.sep+'{}_{}_avg_combined.csv'.format(prename,hm)
            print(cancertype,prename,flag,os.path.isfile(infile),os.path.basename(infile))
            with open(infile) as inf:
                binding_df = pd.read_csv(inf,index_col=0)
            if binding_df.shape[0]>2:
                tmp_df = binding_df.loc[df.index]
                tmp_df.to_csv(outdir+os.sep+'{}_{}_{}.csv'.format(flag,prename,hm))
            #print(tmp_df);exit()
    


def main():

    outdir = 'f1_gained_lost_HM_signals'
    os.makedirs(outdir,exist_ok=True)
    
    cancertype_signamName_match_list = {'T-ALL':['CD4','JURKAT','CUTLL1','PD9','PD31'],'BRCA':['Breast_normal','BRCA'],'LUAD':['Lung_normal','LUAD',],'CRC':['Colon_normal','CRC'],'AML':['Erythroid','AML']}
    cancertype_signamName_match_list = {'AML':['Erythroid','AML']}
    
    for cancertype in cancertype_signamName_match_list.keys():
        gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype)
        abstract_binding_pattern_df(gained_df,cancertype,cancertype_signamName_match_list[cancertype],outdir,'gained')
        abstract_binding_pattern_df(lost_df,cancertype,cancertype_signamName_match_list[cancertype],outdir,'lost')




    
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

