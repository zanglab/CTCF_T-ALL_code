import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import association_with_genes
import association_with_regions
import re,bisect
import CTCF_TALL_modules_new



def abstract_binding_pattern_df(df,cancertype,prenames,outdir,flag):

    combined_file_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f8_panCancer_CTCF_binding_pattern/f2_CTCF_pattern_combined/f1_combined_cancerType_CTCF_pattern_csv'
    for prename in prenames:
        infile = combined_file_dir+os.sep+'{}_avg_combined.csv'.format(prename)
        print(cancertype,os.path.basename(infile))
        with open(infile) as inf:
            binding_df = pd.read_csv(inf,index_col=0)
        tmp_df = binding_df.loc[df.index]
        tmp_df.to_csv(outdir+os.sep+'{}_{}.csv'.format(flag,prename))
        #print(tmp_df);exit()
    


def main():

    outdir = 'f1_gained_lost_CTCF_binding'
    os.makedirs(outdir,exist_ok=True)
    cancertype_signamName_match_list = {'T-ALL':['T-ALL_cancer','T-ALL_normal',\
    'T-ALL_cancer_CUTLL1','T-ALL_cancer_JURKAT','PD9','PD31'],\
    'AML':['AML_cancer','AML_normal'],\
    'BRCA':['BRCA_cancer','BRCA_normal'],\
    'CRC':['CRC_cancer','CRC_normal'],\
    'LUAD':['LUAD_cancer','LUAD_normal'],\
    'PRAD':['PRAD_cancer','PRAD_normal'],\
    'PRAD_TissueAdded':['PRAD_TissueAdded_cancer','PRAD_TissueAdded_normal']}
    
    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']
    for cancertype in cancertypes:
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

