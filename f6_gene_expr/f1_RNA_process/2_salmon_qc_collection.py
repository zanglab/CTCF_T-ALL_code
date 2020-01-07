import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd


def main(indir,outfile):

    salmon_dirs = glob.glob(indir+'/*')
   
    outf = open(outfile,'w')
    outf.write('{}\t{}\t{}\t{}\n'.format('sample_id','total_reads','mapped_reads','mapping_rate'))
    for salmon_dir in sorted(salmon_dirs):
        sample_id = os.path.basename(salmon_dir)
        log_file = salmon_dir+os.sep+'logs'+os.sep+'salmon_quant.log'
        #print(sample_id,os.path.isfile(log_file))
        with open(log_file) as log_f:
            lines = log_f.readlines()
            for line in lines:
                if re.search("Mapping",line):
                    mapping_rate = re.split("\s",line)[-2]
                if re.search("Observed",line):
                    total_reads = re.split("\s",line)[1]      
                if re.search("Counted",line):
                    mapped_reads = re.split("\s",line)[5]                        
            #print(total_reads,mapped_reads,mapping_rate)
            #exit(0)
        outf.write('{}\t{}\t{}\t{}\n'.format(sample_id,total_reads,mapped_reads,mapping_rate))
    outf.close()




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
    parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of salmon results', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<5:
        parser.print_help()
        sys.exit(1)
  
    main(args.indir,args.outfile)
