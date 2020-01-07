import os,sys,argparse
import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
from GenomeData import *
from operator import itemgetter
#def expand_region(summitlist):
import read_write_file 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
# import seaborn as sns
# sns.set(font_scale=1.5)
# sns.set_style("whitegrid", {'axes.grid' : False})


def read_region_from_bed(infile,species,expand,write_out=None,end_col=2,region_sorted=True):
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
                start = int(sline[1]) + int(sline[-1])
                end = start#int(sline[end_col])
                fold_enrichment = float(sline[6])
                if chrom in chroms and fold_enrichment>=4:
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




                       

def main(infile):               

    species='hg38'
    regions = read_region_from_bed(infile,species,expand=0,end_col=1,region_sorted=True)  
    intervals=[] 
    for ele in regions:       
        # calculate the intervals in each chromosome
        mylist = [region[0] for region in regions[ele]]        
        intervals.extend([ mylist[i]-mylist[i-1] for i in np.arange(1,len(mylist))])
                
    #print(len(intervals),min(intervals),max(intervals))
    with open('f1_peak2000_datasets_union_summits/all_summits_fe4_intervals.txt','w') as outf:
        outf.write('\n'.join([str(i) for i in intervals]))


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of bed fromat, union all the overlapping regions', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)

    main(args.infile)