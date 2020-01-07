import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
# import association_with_genes
# import association_with_regions
import re,bisect
# import CTCF_TALL_modules
import scipy
from scipy import stats
#import find_overlap_keep_info_NOT_sep_strand_lastColMarked

tf_455_file = "/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2017/tf_merged/log_tf_binding_total_hg38.txt"
# tf_455_file = "/Volumes/zanglab/zw5j/data/cistrome/cistrome_db_2017/tf_merged/log_tf_binding_total_hg38.txt"
#tf_455_file="/nv/vol190/zanglab/zw5j/data/ChIPseq_overlap_2016_new/analysis/Info_455/tf_455_total.txt"
tf_455_df = pd.read_csv(tf_455_file,sep='\t',index_col=0,header=None)
tf_455_df.columns = ['udhs_binding']
udhs_total = 2723010

def stats_test(treat_df,basename,outdir):

    df = pd.DataFrame()
    tf_columns = treat_df.columns[1:]
    treat_total = treat_df.shape[0]
    for motif in tf_columns:
        treat_nonzero = len([i for i in treat_df[motif] if i !=0])
        udhs_nonzero = tf_455_df.loc[motif,'udhs_binding']
        s3,p3 = stats.fisher_exact([[treat_nonzero,treat_total-treat_nonzero],[udhs_nonzero,udhs_total-udhs_nonzero]])
        df.loc[motif,'t_binding'] = treat_nonzero
        df.loc[motif,'t_total'] = treat_total
        df.loc[motif,'c_binding'] = udhs_nonzero
        df.loc[motif,'c_total'] = udhs_total
        df.loc[motif,'udhs_ctrl_fisher_s'] = np.round(s3,2)
        df.loc[motif,'udhs_ctrl_fisher_p'] = '{:.2e}'.format(p3)

    df = df.sort_values(by=['udhs_ctrl_fisher_s'],ascending=False)
    df.to_csv(outdir+os.sep+'{}_fisher.csv'.format(basename))
    

def stats_test_with_control(treat_df,ctrl_df,basename,outdir):

    df = pd.DataFrame()
    
    tf_columns = treat_df.columns[1:]  
    treat_total = treat_df.shape[0]
    ctrl_total = ctrl_df.shape[0]
    for motif in tf_columns:
        treat_nonzero = len([i for i in treat_df[motif] if i !=0])
        ctrl_nonzero = len([i for i in ctrl_df[motif] if i !=0])
        s,p = stats.fisher_exact([[treat_nonzero,treat_total-treat_nonzero],[ctrl_nonzero,ctrl_total-ctrl_nonzero]])
        df.loc[motif,'t_binding'] = treat_nonzero
        df.loc[motif,'t_total'] = treat_total
        df.loc[motif,'c_binding'] = ctrl_nonzero
        df.loc[motif,'c_total'] = ctrl_total
        df.loc[motif,'ctrl_fisher_s'] = np.round(s,2)
        df.loc[motif,'ctrl_fisher_p'] = '{:.2e}'.format(p)

    df = df.sort_values(by=['ctrl_fisher_s'],ascending=False)
    df.to_csv(outdir+os.sep+'{}_with_ctrl_fisher.csv'.format(basename))



def main(treat,ctrl,outdir):

    basename=os.path.basename(treat).split('_udhs_chip_peak_overlapped.csv')[0];print(basename)
    treat_df = pd.read_csv(treat,index_col=0)
    control_df = pd.read_csv(ctrl,index_col=0)
    #print(control_df)
    stats_test_with_control(treat_df,control_df,basename,outdir)

        

 
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--treat', action = 'store', type = str,dest = 'treat', help = 'input file of', metavar = '<file>')
    parser.add_argument('-c','--control', action = 'store', type = str,dest = 'control', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.treat,args.control,args.outdir)

