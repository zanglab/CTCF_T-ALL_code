import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
# import multiprocessing


def get_lines(infile):

    with open(infile,'rb') as f:
        lines = 0
        buf_size = 1024*1024
        buf = f.raw.read(buf_size)
        while buf:
            lines += buf.count(b'\n')
            buf = f.raw.read(buf_size)
    return lines



def main(infile):

    df = pd.read_csv('../f1_GSM2SRR/ctcf_new_GSM_to_SRR.csv',sep=',',index_col=0)    
    peak_count_df = pd.DataFrame()
    for gsmid in df.index:
        peak_file = 'GSM_callpeak_out/{}/{}_peaks.narrowPeak'.format(gsmid,gsmid)
        if os.path.isfile(peak_file):
            num_peaks = get_lines(peak_file)
            peak_count_df.loc[gsmid,'#peaks'] = num_peaks
        else:
            print('rm GSM_callpeak_out/{}/{}*'.format(gsmid,gsmid))
#             os.system('rm GSM_callpeak_out/{}/{}*'.format(gsmid,gsmid))#;exit()
    peak_count_df.to_csv("GSM_CTCF_peaks.csv")
        



    
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
  
    main(args.infile)
