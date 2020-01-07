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


cancer_type_matched = {'T-ALL':['CD4','JURKAT','CUTLL1','PD9','PD31'],'BRCA':['Breast_normal','BRCA'],'LUAD':['Lung_normal','LUAD',],'CRC':['Colon_normal','CRC'],'PRAD':['CRC'],'AML':['Erythroid','AML']}
# cancer_type_matched = {'AML':['Erythroid']}
   
def main():

    outdir = 'f1_combined_cancerType_csv'
    os.makedirs(outdir,exist_ok=True)

    for cancertype in cancer_type_matched.keys():
        for celltype in cancer_type_matched[cancertype]:
            prename_df = pd.read_excel('/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f9_histone_modification/data_file/CTCF_pancancer_histone_modification_201907.xlsx',sheet_name=celltype)
            #print(celltype,df)#;exit()
            for hm in ['H3K27ac','H3K27me3','H3K4me1']:
                prenames = prename_df[hm].dropna()
                df = pd.DataFrame()
                i=0
                for prename in prenames:
                    try:
                        prename=int(prename)
                    except:
                        pass
                    csv_file = '{}/{}_{}_{}_{}.txt'.format('../f1_binding_pattern_on_union/f1_HM_pattern_csv/',cancertype,celltype,hm,prename)
                    if os.path.isfile(csv_file):
                        with open(csv_file) as csv_inf:
                            sample_df = pd.read_csv(csv_inf,sep="\t",index_col=0)#;print(sample_df);exit()
                        if i==0:
                            df = sample_df
                        else:
                            df = df+sample_df
                        i+=1
                    #print(sample_df,df)
                df = df/i 
                df = df.round(2)     
                #print(df);exit()

                df.to_csv(outdir+os.sep+'{}_{}_avg_combined.csv'.format(celltype,hm))
 

    



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
