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

    outdir = 'f2_combined_HM_RPKM_csv'
    os.makedirs(outdir,exist_ok=True)

    celltypes=['CD4','JURKAT']
    for celltype in celltypes:
        prename_df = pd.read_excel('/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f11_histone_modification/data_files/CTCF_pancancer_histone_modification.xlsx',sheet_name=celltype)
            
        for hm in ['H3K27ac']:
            prenames = prename_df[hm].dropna()
            df = pd.DataFrame()
            i=0
            for prename in prenames:
                try:
                    prename=int(prename)
                except:
                    pass
                csv_file = '{}/{}_{}_{}.txt'.format('f1_HM_RPKM/',celltype,hm,prename)
                #print(os.path.isfile(csv_file))

                with open(csv_file) as csv_inf:
                    sample_df = pd.read_csv(csv_inf,sep="\t",index_col=0,header=None)#;print(sample_df);exit()
                if i==0:
                    df = sample_df
                else:
                    df = pd.concat([df,sample_df],axis=1)
                i+=1
                #if i==2:
                 #   print(sample_df,df);exit()
            #df = df/i 
            df = df.round(2)     
            #print(df);exit()

            df.to_csv(outdir+os.sep+'{}_{}_combined.csv'.format(celltype,hm))
 

    



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
