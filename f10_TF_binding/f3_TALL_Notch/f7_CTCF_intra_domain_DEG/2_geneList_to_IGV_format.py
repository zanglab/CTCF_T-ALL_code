import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import CTCF_TALL_modules_new
import scipy
from scipy import stats
import association_with_genes   


def generate_IGV_format(gene_list,ucsc_df,outdir):
    basename = os.path.basename(gene_list).split('.txt')[0]
    genes = [i.strip() for i in open(gene_list).readlines()]
    outf = open(outdir+os.sep+basename+'_IGV.bed','w')
    for gene in ucsc_df.index.intersection(genes):
        chr = ucsc_df.loc[gene].chrom
        strand = ucsc_df.loc[gene].strand
        gbstart = ucsc_df.loc[gene].gbStart
        gbend = ucsc_df.loc[gene].gbEnd
        outf.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr,gbstart,gbend,gene,'0',strand))
    outf.close()




def main():

    outdir = 'f2_genes_IGV'
    os.makedirs(outdir,exist_ok=True) 

    ucsc_df = association_with_genes.return_ucsc_df('hg38')
    gene_lists = glob.glob('f1_intra_domain_deg/*txt')
    for gene_list in gene_lists:
        generate_IGV_format(gene_list,ucsc_df,outdir)




 
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

