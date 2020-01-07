import os,sys,argparse,glob
import numpy as np
import pandas as pd
from GenomeData import *
import re,bisect
from operator import itemgetter

plus = re.compile("\+")
minus = re.compile("\-")

def union_islands(islandlist):
    #start/[0] position in islandlist must be sorted
    unionlist = []
    i = 1
    currentStart,currentEnd = islandlist[0][0],islandlist[0][1]
    while i<len(islandlist):
        compareStart,compareEnd = islandlist[i][0],islandlist[i][1]
        assert currentStart <= compareStart
        if compareStart > currentEnd:
            unionlist.append([currentStart,currentEnd])
            currentStart,currentEnd = compareStart,compareEnd 
            i+=1
        else:
            currentEnd = max(currentEnd,compareEnd)
            i+=1    
        #print(i)   
    unionlist.append([currentStart,currentEnd])
    return unionlist
    
def get_overlap_info_strID(islands,regions):
    # for each island in islands[chrom], check if 1-overlap 0-non-overlap with regions[chrom]
    # islands/regions represents for compare_regions/basic_region respectively
    #island_dict = {}
    #i=0
    overlapped,nonoverlapped = {},{}
    for chrom in islands:
        if chrom not in regions:
            for island in islands[chrom]:
                if chrom not in nonoverlapped:
                    nonoverlapped[chrom] = []
                nonoverlapped[chrom].append(island)
        else:
            # check the regions in regions[chrom] are non-overlap and sorted
            regionlist = []
            for region in regions[chrom]:
                regionlist.extend(region)
            for num in range(len(regionlist)-1):
                assert regionlist[num]<regionlist[num+1]
                            
            # for each island, check if overlapped with regions[chrom]
            for island in islands[chrom]:
                #ID = chrom+'_'+str(island[0])+'_'+str(island[1])+'_'+str(i)
                assert island[0]<=island[1]
                s = bisect.bisect_left(regionlist,island[0])
                e = bisect.bisect_right(regionlist,island[1])
                if s==e and s%2==0:
                    if chrom not in nonoverlapped:
                        nonoverlapped[chrom] = []
                    nonoverlapped[chrom].append(island)
                else:
                    if chrom not in overlapped:
                        overlapped[chrom] = []
                    overlapped[chrom].append(island)
                
    return overlapped,nonoverlapped

def get_overlap_info_numID(islands,regions):
    # for each island in islands[chrom], check if 1-overlap 0-non-overlap with regions[chrom]
    overlapped,nonoverlapped = {},{}
    for chrom in islands:
        
        if chrom not in regions:
            for island in islands[chrom]:
                if chrom not in nonoverlapped:
                    nonoverlapped[chrom] = []
                nonoverlapped[chrom].append(island)
        else:
            # check the regions in regions[chrom] are non-overlap and sorted
            regionStart = [region[0] for region in regions[chrom]]
            regionEnd = [region[1] for region in regions[chrom]]
                   
            # for each island, check if overlapped with regions[chrom]
            for island in islands[chrom]:              
                assert island[0]<=island[1]
                e = bisect.bisect_left(regionEnd,island[0])
                if e==0:
                    left_pos = 0
                    right_pos = regions[chrom][e][0]
                elif e==len(regionEnd):
                    left_pos = regions[chrom][e-1][0]
                    right_pos = 999999999
                else:
                    left_pos = regions[chrom][e-1][0]
                    right_pos = regions[chrom][e][0]

                island[-1]=island[-1]+'\t{}\t{}'.format(left_pos,right_pos)
                if chrom not in overlapped:
                    overlapped[chrom] = []
                overlapped[chrom].append(island)
                #i+=1
    return overlapped,nonoverlapped

def read_region_from_bed(infile,species,expand=0,write_out=None,end_col=2,region_sorted=True):
    if species in species_chroms.keys():
        chroms = species_chroms[species]
    else:
        print('Species not recognized!')
        exit(1)
    regions1,regions2 = {},{}
    with open(infile,'r') as inf:
        line = inf.readline()
        while line:
            #print(line)
            if not re.match("#",line):
                #print(line)
                sline = line.strip().split('\t')
                chrom = sline[0]
                start = int(sline[1])
                end = int(sline[end_col])
                #strand = sline[5]
                #if plus.match(strand):
                if 1:
                    if chrom in chroms:
                        if chrom not in regions1:
                            regions1[chrom]=[]
                        #regions1[chrom].append([max(0,start-expand),end+expand,'\t'.join(sline[3:])])
                        mid=int((start+end)*.5)
                        regions1[chrom].append([mid,mid+1,'\t'.join(sline[3:])])
                    else:
                        pass   
                elif minus.match(strand):
                    if chrom in chroms:
                        if chrom not in regions2:
                            regions2[chrom]=[]
                        regions2[chrom].append([max(0,start-expand),end+expand,'\t'.join(sline[3:])])
                    else:
                        pass
                                                       
            #print(regions1,regions2);exit(0)
            line = inf.readline()
            
    if region_sorted:
        for chrom in regions1:
            regions1[chrom] = sorted(regions1[chrom],key = itemgetter(0),reverse=False)
        for chrom in regions2:
            regions2[chrom] = sorted(regions2[chrom],key = itemgetter(0),reverse=False)
            
    if write_out:
        with open(write_out,'w') as outf:
            for chrom in regions1:
                for region in regions1[chrom]:
                    outf.write('{}\t{}\t{}\t{}\n'.format(chrom,region[0],region[1],region[2]))
            for chrom in regions2:
                for region in regions2[chrom]:
                    outf.write('{}\t{}\t{}\t{}\n'.format(chrom,region[0],region[1],region[2]))
                        
    return regions1,regions2


def write_regions_bed3(regions1,regions2,outfile):
    # write the regions(with chrom,start,end info) into bed3 format
    with open(outfile,'w') as outf:
        outf.write('{}\t{}\t{}\t{}\n'.format('#chrom','start','end','id\tscore\tstrand\tleft_pos\tright_pos'))
        for chrom in regions1:
            for region in regions1[chrom]:
                outf.write('{}\t{}\t{}\t{}\n'.format(chrom,region[0],region[1],region[2]))
        for chrom in regions2:
            for region in regions2[chrom]:
                outf.write('{}\t{}\t{}\t{}\n'.format(chrom,region[0],region[1],region[2]))  


def main(infile1,infile2,species,outfile1,expand1,expand2):

    compare_regions1,compare_regions2 = read_region_from_bed(infile1,species,expand1)
    basic_regions1,basic_regions2 = read_region_from_bed(infile2,species,expand2)
    basic_regions_union = {}
    for ele in basic_regions1:
        basic_regions_union[ele] = union_islands(basic_regions1[ele])
    #for ele in basic_regions2:
        #basic_regions2[ele] = union_islands(basic_regions2[ele])
    
    # test if #union regions equals to #non union regions
    a = sum(len(basic_regions1[i]) for i in basic_regions1)
    b = sum(len(basic_regions_union[i]) for i in basic_regions_union)
    assert a==b
    
    overlapped1,nonoverlapped1 = get_overlap_info_numID(compare_regions1,basic_regions1)        
    #overlapped2,nonoverlapped2 = get_overlap_info_strID(compare_regions2,basic_regions2)
    print(nonoverlapped1)
    #if outfile1:
    write_regions_bed3(overlapped1,nonoverlapped1,outfile1) 
    #if outfile2:
    #write_regions_bed3(nonoverlapped1,nonoverlapped2,outfile2)  
    
    #test_df = pd.read_csv('dynamic_notch_co_binding/dynamic_notch_domain_co_occurrence.bed',sep='\t',index_col=4) 
    #print(test_df)      
        
        
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--infile1', action = 'store', type = str,dest = 'infile1', help = 'input file to be compared/separated', metavar = '<file>')
    parser.add_argument('-b', '--infile2', action = 'store', type = str,dest = 'infile2', help = 'input file to be compared as basic', metavar = '<file>')
    parser.add_argument('-p','--outfile1', action = 'store', type = str,dest = 'outfile1', help = 'outfile of 1-overlap/0-nonoverlap info with b', metavar = '<file>')
    #parser.add_argument('-q','--outfile2', action = 'store', type = str,dest = 'outfile2', help = 'outfile of a not overlapped with b', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    parser.add_argument('-e1', '--expand1', action = 'store', type = int,dest = 'expand1', help = 'expand of file1', metavar = '<file>',default=0)
    parser.add_argument('-e2', '--expand2', action = 'store', type = int,dest = 'expand2', help = 'expand of file2', metavar = '<file>',default=0)
    

    args = parser.parse_args()
    if(len(sys.argv))<9:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile1,args.infile2,args.species,args.outfile1,args.expand1,args.expand2)
