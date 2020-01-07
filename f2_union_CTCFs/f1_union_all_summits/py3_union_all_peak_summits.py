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
import association_with_regions


def read_summits_from_narrowpeak_bed(infile,species,expand,write_out=None,end_col=2,region_sorted=True):
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

    


def write_regions_bed3(regions,outfile):
    # write the regions(with chrom,start,end info) into bed3 format
    with open(outfile,'w') as outf:
        i = 0
        for chrom in regions:
            for region in regions[chrom]:
                outf.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom,region[0],region[1],i,0,'+'))
                i = i+1

def main():               
    
    peak_file='f1_peak2000_datasets_union_summits/collection_all_CTCF_narrowPeak.bed'
    species = 'hg38'
    # union all peaks, 150bp expanded summits
    for width in [150]:
        regions = read_summits_from_narrowpeak_bed(peak_file,species,expand=int(width/2),end_col=1,region_sorted=True)   
        for ele in regions:
            regions[ele] = association_with_regions.union_islands(regions[ele])       
        write_regions_bed3(regions,'f1_peak2000_datasets_union_summits/union_summits_fe4_width_{}.bed'.format(width))
        


            
            
if __name__=='__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of bed fromat, union all the overlapping regions', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)

    main()