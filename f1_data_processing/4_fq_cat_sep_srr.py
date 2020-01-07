import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
# import multiprocessing





def main(infile):

    indir='../f3_fq_files'
     
    df = pd.read_csv('f1_GSM2SRR/ctcf_new_GSM_to_SRR.csv',sep=',',index_col=0)    
    for gsmid in df.index:
#         print('##',gsmid)
        srrs = df.loc[gsmid,'SRR'].split('&&')
        layout = df.loc[gsmid,'layout']
        if len(srrs)!=1:
            outf = open('run_py4_slurms/{}.slurm'.format(gsmid),'w')
            outf.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
''')
            outf.write('#SBATCH -o out_{}.log'.format(gsmid))
            outf.write('\n\n##{}'.format(gsmid))
            if layout=='SINGLE':
                # == cat *.fastq files ==
                outf.write('\n\ncat ')
                for SRRid in srrs:
                    fq_file='{}/{}.fastq'.format(indir,SRRid)
                    outf.write('{} '.format(fq_file))
                outf.write(' > ../f3_fq_files/{}.fastq\n'.format(gsmid))
                outf.write('wait\ngzip ../f3_fq_files/{}.fastq'.format(gsmid))
            else:
                # == cat *_1.fastq files ==
                outf.write('\n\ncat ')
                for SRRid in srrs:
                    fq_file1='{}/{}_1.fastq'.format(indir,SRRid)
                    outf.write('{} '.format(fq_file1))
                outf.write('> ../f3_fq_files/{}_1.fastq\n'.format(gsmid))
                outf.write('wait\ngzip ../f3_fq_files/{}_1.fastq'.format(gsmid))
                # == cat *_2.fastq files ==
                outf.write('\n\ncat ')
                for SRRid in srrs:
                    fq_file2='{}/{}_2.fastq'.format(indir,SRRid)
                    outf.write('{} '.format(fq_file2))
                outf.write('> ../f3_fq_files/{}_2.fastq\n'.format(gsmid))
                outf.write('wait\ngzip ../f3_fq_files/{}_2.fastq'.format(gsmid))
#             exit()

#     print('\n\n#### rename files ####')
#     for gsmid in df.index:
#         srrs = df.loc[gsmid,'SRR'].split('&&')
#         layout = df.loc[gsmid,'layout']
#         if len(srrs)==1:
#             print('\n##',gsmid)
#             if layout=='SINGLE':
#                 # == rename *.fastq files ==
#                 print('mv',end=' ')
#                 for SRRid in srrs:
#                     fq_file='{}/{}.fastq.gz'.format(indir,SRRid)
#                     print(fq_file,end=' ')
#                 print('f3_fq_files/{}.fastq.gz'.format(gsmid))
#             else:
#                 # == rename *_1.fastq files ==
#                 print('mv',end=' ')
#                 for SRRid in srrs:
#                     fq_file1='{}/{}_1.fastq.gz'.format(indir,SRRid)
#                     print(fq_file1,end=' ')
#                 print('f3_fq_files/{}_1.fastq.gz'.format(gsmid))
#                 
#                 # == rename *_2.fastq files ==
#                 print('mv',end=' ')
#                 for SRRid in srrs:
#                     fq_file2='{}/{}_2.fastq.gz'.format(indir,SRRid)
#                     print(fq_file2,end=' ')
#                 print('f3_fq_files/{}_2.fastq.gz'.format(gsmid))



    
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
