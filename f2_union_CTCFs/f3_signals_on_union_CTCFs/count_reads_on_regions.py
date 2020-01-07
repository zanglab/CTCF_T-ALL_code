# Time-stamp: <2017-08-10>
'''
Copyright (c) 2017, 2018 Chongzhi Zang, Zhenjia Wang <zhenjia@virginia.edu>

This code is free software; you can redistribute it and/or modify it 
under the terms of the BSD License.

@status: release candidate
@version: $Id$
@author: Chongzhi Zang, Zhenjia Wang
@contact: zhenjia@virginia.edu

This file is used to count the reads mapped onto bed format regions
'''

import os,re,argparse,sys
import bisect
from IOparser_BedBam import get_tag_regions
#from BART.OptValidator import

plus = re.compile('\+')
minus = re.compile('\-')
      

def is_list_sorted(mylist):
    '''
    Check if list is sorted
    '''        
    for i in range(len(mylist)-1):
        if mylist[i] > mylist[i+1]:
            return 0
    return 1   
        	

def get_read_positions(positions,regions,val,fragment_size):
    '''
    Return all the shifted positions of reads 
    '''
    # default fragment_size = 150
    if fragment_size >= 0:
        shift = int(round(fragment_size/2))	
        for chrom in regions.keys():
            if chrom not in positions:
                positions[chrom]=[]
            for inner in regions[chrom].keys():
                positions[chrom].extend([outer+shift*val for outer in regions[chrom][inner]])
    else:
        for chrom in regions.keys():
            if chrom not in positions:
                positions[chrom]=[]
            for inner in regions[chrom].keys():
                for outer in regions[chrom][inner]:
                    shift = int(abs(outer-inner)/2)	
                    positions[chrom].append(outer+shift*val)
    return positions   
	        	    
		
def get_count_on_mapID(start,end,positions,mid,expand):
    '''
    Count the tags/positions on DHS
    '''
    #if is_list_sorted(positions)==0:
    #    positions.sort()
    if start < end and mid:
        middle = int((start+end)*0.5)
        s = bisect.bisect_left(positions,middle-expand)
        e = bisect.bisect_right(positions,middle+expand)
        return e-s, expand*2
        
    elif start < end:
        s = bisect.bisect_left(positions,start-expand)
        e = bisect.bisect_right(positions,end+expand)
        return e-s, end-start+expand*2
            
    else:
        return 0,1


def read_count_on_mapfile(infile,readsfile,species,format,fragmentsize=147,mid=False,expand=0):
    '''
    Count the num of (unique) reads in user-input bed/bam file
    on each of the UDHS -- bart profile
    '''
    #specify the species
    # get the start-end regions of each read in each chrom (as key) and separate by strand

    regions1, regions2 = get_tag_regions(species,format,readsfile)

    # get the counting position of each read(tag), chrom as key
    positions = get_read_positions({},regions1,1,fragmentsize)
    positions = get_read_positions(positions,regions2,-1,fragmentsize)
    total=0
    for chrom in positions:
        positions[chrom].sort() 
        total+=len(positions[chrom])
    # check if positions is NULL
    #if total ==0:
        #sys.stderr.write('Can not read the input bed/bam file!\n')
    # count reads on each DHS

    
    mapfile = open(infile,'r')
    line = mapfile.readline()
    counting = {}
    while line:
        line = line.strip().split()
        #dhs = BED(line[0],line[1],line[2],line[3],line[4],line[5])
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        map_id = line[3]
        if chrom in positions:
            #assert DHS_id not in counting
            nums,region_width = get_count_on_mapID(start,end,positions[chrom],mid,expand)
            counting[map_id]= round(nums*1000000000/(total*region_width),3)
        line = mapfile.readline()
    mapfile.close()
    return counting

def main(infile,readsfile,outfile,species,format,fragmentsize,mid,expand):

    counting = read_count_on_mapfile(infile,readsfile,species,format,fragmentsize,mid,expand)
    with open(outfile,'w') as outf:
        for i in counting:
            outf.write('{}\t{}\n'.format(i,counting[i]))
    
    
    
    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    parser.add_argument('-b', '--readsfile', action = 'store', type = str,dest = 'readsfile', help = 'ChIP-seq reads in bed/bam format', metavar = '<file>')
    parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile to write the the RPKM for each region', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
   
    parser.add_argument('-f', '--format', action = 'store', type = str,dest = 'format', help = 'format of file, BED or BAM', metavar = '<str>',required=True)
    parser.add_argument('-g', '--fragmentsize', action = 'store', type = int,dest = 'fragmentsize', help = 'fragmentsize for the shift of reads. Default: 147', metavar = '<int>',default=147)
    parser.add_argument('-m', '--mid', action = 'store_true', dest = 'mid', help = 'whether to use middle side for expansion. Default: False',default=False)
    parser.add_argument('-e', '--expand', action = 'store', type = int,dest = 'expand', help = 'expand of regions. Default: 0', metavar = '<int>',default=0)
    

    args = parser.parse_args()
    if(len(sys.argv))<9:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile,args.readsfile,args.outfile,args.species,args.format,args.fragmentsize,args.mid,args.expand)
