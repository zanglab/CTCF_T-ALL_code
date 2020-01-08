import sys,argparse
from GenomeData import *
from operator import itemgetter
#import os,glob
#import numpy as np
#import pandas as pd
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

#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')


def get_lines(infile):

    with open(infile,'rb') as f:
        lines = 0
        buf_size = 1024*1024
        buf = f.raw.read(buf_size)
        while buf:
            lines += buf.count(b'\n')
            buf = f.raw.read(buf_size)
    return lines


def get_list_of_lines(infile,*args):

    return len(args)

def read_region_from_bed(infile,species,expand=0,write_out=None,end_col=2,region_sorted=True):
    if species in species_chroms.keys():
        chroms = species_chroms[species]
    else:
        print('Species not recognized!')
        exit(1)
    regions = {}
    with open(infile,'r') as inf:
        line = inf.readline()
        while line:
            #print(line)
            if not re.match("#",line):
                #print(line)
                sline = line.strip().split()
                chrom = sline[0]
                start = int(sline[1])
                end = int(sline[end_col])
                if chrom in chroms:
                    if chrom not in regions:
                        regions[chrom]=[]
                    regions[chrom].append([max(0,start-expand),end+expand])
                else:
                    pass                   
            line = inf.readline()
    if region_sorted:
        for chrom in regions:
            regions[chrom] = sorted(regions[chrom],key = itemgetter(0),reverse=False)
    if write_out:
        with open(write_out,'w') as outf:
            for chrom in regions:
                for region in regions[chrom]:
                    outf.write('{}\t{}\t{}\n'.format(chrom,region[0],region[1]))
    
    return regions
                
def write_regions_bed3(regions,outfile):
    # write the regions(with chrom,start,end info) into bed3 format
    with open(outfile,'w') as outf:
        for chrom in regions:
            for region in regions[chrom]:
                outf.write('{}\t{}\t{}\n'.format(chrom,region[0],region[1]))
    
    
 
def write_counter_table(counter_table,outfilename,plot=False):

    with open(outfilename+'.txt','w') as outf:
        for counter_key in sorted(counter_table.keys(),key=int):
            outf.write('{}\t{}\n'.format(counter_key,counter_table[counter_key]))
    if plot:
        fig = plt.figure()
        x = sorted(counter_table.keys(),key=int)
        y = [counter_table[counter_key] for counter_key in sorted(counter_table.keys(),key=int)]
        plt.bar(x,y)
        plt.savefig(outfilename+'.png',bbox_inches='tight')
        plt.close()
        

def main(infile):
    pass

    




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
#    parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
#    parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<3:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile,args.species)
