'''
Created on 10/04/2018
@authors: Chongzhi Zang, Zhenjia Wang<zw5j@virginia.edu>

This file is used to get methylation state on each region/position
-- centered by middle with given expand value, or
-- start to end position

'''

import os,re,argparse,sys
import bisect
from IOparser_BedBam import get_tag_regions
#from BART.OptValidator import
import association_with_regions
import numpy as np
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
        	
		
def get_count_on_mapID(chrom,start,end,tags_startlist,tags_regions,mid,expand,colnums):
    '''
    Count the tags/positions on DHS
    '''
    #if is_list_sorted(positions)==0:
    #    positions.sort()
    
    states = [-1]*(end-start+1)
    s = bisect.bisect_left(tags_startlist,start)
    e = bisect.bisect_right(tags_startlist,end)

    if s!=e:
        for i in np.arange(s,e):
            score_list = [float(val) for val in tags_regions[i][-1].split('\t')[-1*colnums:]]
            #print(start,end,tags_regions[i],score_list);exit()
            pos = tags_regions[i][0]-start
            states[pos] = np.mean(score_list)

    return states


def read_count_on_mapfile(infile,readsfile,outfile,species,mid=False,expand=0,colnums=1):
    '''
    Count the num of (unique) reads in user-input bed/bam file
    on each of the UDHS -- bart profile
    '''
    #specify the species
    # get the start-end regions of each read in each chrom (as key) and separate by strand

    tags_regions = association_with_regions.read_regions_from_bed_not_sep_strand(readsfile,species)
    tags_startlist = {}
    for ele in tags_regions:
        tags_startlist[ele] = [i[0] for i in tags_regions[ele]]
    
    outf = open(outfile,'w')
    mapfile = open(infile,'r')
    line = mapfile.readline()
    while line:
        line = line.strip().split()
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        map_id = line[3]
        if start < end and mid:
            middle = int((start+end)*0.5)
            start = middle-expand
            end = middle+expand
        elif start < end:
            start = start-expand
            end = end+expand
        
        states = [-1]*(end-start+1)
        if chrom in tags_startlist:
            states = get_count_on_mapID(chrom,start,end,tags_startlist[chrom],tags_regions[chrom],mid,expand,colnums)
        outf.write('{}\t{}\n'.format(map_id,','.join([str(i) for i in states])))
        line = mapfile.readline()
    mapfile.close()
    outf.close()


def main(infile,readsfile,outfile,species,mid,expand,colnums):

    read_count_on_mapfile(infile,readsfile,outfile,species,mid,expand,colnums)
    #with open(outfile,'w') as outf:
    #    for i in counting:
    #        outf.write('{}\t{}\n'.format(i,','.join([str(i) for i in counting[i]])))
    
    
    
    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    parser.add_argument('-t', '--tagsfile', action = 'store', type = str,dest = 'tagsfile', help = 'ChIP-seq reads in bed/bam format', metavar = '<file>')
    parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile to write the the RPKM for each region', metavar = '<file>',required=True)
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
   
    #parser.add_argument('-f', '--format', action = 'store', type = str,dest = 'format', help = 'format of file, BED or BAM', metavar = '<str>',required=True)
    #parser.add_argument('-g', '--fragmentsize', action = 'store', type = int,dest = 'fragmentsize', help = 'fragmentsize for the shift of reads. Default: 147', metavar = '<int>',default=147)
    parser.add_argument('-m', '--mid', action = 'store_true', dest = 'mid', help = 'whether to use middle side for expansion. Default: False',default=False)
    parser.add_argument('-e', '--expand', action = 'store', type = int,dest = 'expand', help = 'expand of regions. Default: 0', metavar = '<int>',default=0)
    parser.add_argument('-c', '--colnums', action = 'store', type = int,dest = 'colnums', help = 'num of -1*cols with scores. Default: 1 (0 for all chr-start-end excluede columns, e.g., in union bdg format)', metavar = '<int>',default=1)    
    

    args = parser.parse_args()
    if(len(sys.argv))<9:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile,args.tagsfile,args.outfile,args.species,args.mid,args.expand,args.colnums)
