import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
# import multiprocessing



def main(infile):

     
    df = pd.read_csv('../f1_GSM2SRR/ctcf_new_GSM_to_SRR.csv',sep=',',index_col=0)    
    for gsmid in df.index:
        srrs = df.loc[gsmid,'SRR'].split('&&')
        layout = df.loc[gsmid,'layout']
        if layout=='SINGLE':
            print('++++\n#name:\n{}\n#treat'.format(gsmid))
            fq_file='{}.fastq.gz'.format(gsmid)#;print(os.path.isfile(fq_file))
            print('{}\n'.format(fq_file))
        
    for gsmid in df.index:
        srrs = df.loc[gsmid,'SRR'].split('&&')
        layout = df.loc[gsmid,'layout']
        if layout!='SINGLE':
            print('++++\n#name:\n{}\n#treat'.format(gsmid))
            fq_file1='{}_1.fastq.gz'.format(gsmid)
            fq_file2='{}_2.fastq.gz'.format(gsmid)#;print(os.path.isfile(fq_file2))
            print('{}\n{}\n'.format(fq_file1,fq_file2))
            





    
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
