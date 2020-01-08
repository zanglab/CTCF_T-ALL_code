'''
Created on 09 Dec, 2017

@authors: Zhenjia Wang(zhenjia@virginia.edu), Chongzhi Zang

This file is used to return region/tags associations

modules including:
	union_islands
	fileter_BED_by_regions_not_sep_strand
	fileter_BED_by_regions_not_sep_strand_keep_score
	fileter_BED_by_regions_sep_strand
	get_overlap_info_regions2regions_strID
	get_overlap_info_regions2regions_numID
	get_overlap_info_BEDtag2regions
	read_regions_from_bed_sep_strand
	read_regions_from_bed_not_sep_strand
	read_regions_from_bed_by_chrom_not_sep_strand
	regionslist_from_regions
	...

'''

import time
import os,sys,argparse,glob
import numpy as np
import pandas as pd
from GenomeData import *
import re,bisect
from operator import itemgetter

plus = re.compile("\+")
minus = re.compile("\-")

def sort_islands(islandlist):
    return sorted(islandlist,key = itemgetter(0),reverse=False)
    
def union_islands(islandlist,expand=0):
    #start/[0] position in islandlist must be sorted
    unionlist = []
    i = 1
    try:
        currentStart,currentEnd = islandlist[0][0]-expand,islandlist[0][1]+expand
    except:
        print(islandlist)
    while i<len(islandlist):
        compareStart,compareEnd = islandlist[i][0]-expand,islandlist[i][1]+expand
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

def union_islands_extend_ids(islandlist,expand=0):
    #start/[0] position in islandlist must be sorted
    unionlist = []
    i = 1
    try:
        currentStart,currentEnd,currentId = islandlist[0][0]-expand,islandlist[0][1]+expand,[islandlist[0][2]]
    except:
        print(islandlist)
    while i<len(islandlist):
        compareStart,compareEnd,compareId = islandlist[i][0]-expand,islandlist[i][1]+expand,[islandlist[i][2]]
        assert currentStart <= compareStart
        if compareStart > currentEnd:
            unionlist.append([currentStart,currentEnd,currentId])
            currentStart,currentEnd,currentId = compareStart,compareEnd ,compareId
            i+=1
        else:
            currentEnd = max(currentEnd,compareEnd)
            currentId.extend(compareId)
            i+=1    
        #print(i)   
    unionlist.append([currentStart,currentEnd,currentId])
    return unionlist
    
def union_regions_from_2strands(regionsplus,regionsminus):

    unionregions = regionsplus.copy()
    for chrom in regionsminus:
        if chrom in unionregions:
            unionregions[chrom].extend(regionsminus[chrom])
            #unionregions[chrom] = sorted(unionregions[chrom],key = itemgetter(0),reverse=False)
        else:
            unionregions[chrom] = regionsminus[chrom]
    for chrom in unionregions:
        unionregions[chrom] = sorted(unionregions[chrom],key = itemgetter(0),reverse=False)
    return unionregions

def region_start_end_list_with_expansion(regions,expand = 0):
    region_start_list,region_end_list = {},{}
    for chr in regions:
        region_start_list[chr] = [i[0]-expand for i in regions[chr]]
        region_end_list[chr] = [i[1]+expand for i in regions[chr]]
    return region_start_list,region_end_list

def region_start_end_list(regions):
    region_start_list,region_end_list = {},{}
    for chr in regions:
        region_start_list[chr] = [i[0] for i in regions[chr]]
        region_end_list[chr] = [i[1] for i in regions[chr]]
    return region_start_list,region_end_list
    

def regionslist_from_regions(regions):
    # return the s-e-s-e.... regions list for all start-end regions
    regionlist = []
    for region in regions:
        regionlist.extend(region)
    for num in range(len(regionlist)-1):
        assert regionlist[num]<=regionlist[num+1]  
    return regionlist  

def get_overlap_info_BEDtag2regions(tag,regions,info_column=False):
    # if the file is too big, compare tags one by one
    # each regions[chr] is regionslist from regionslist_from_regions(regions)
    # return tag[column] if necessary
    chrom = tag[0]
    start = int(tag[1])
    end = int(tag[2])
    if chrom not in regions:
        return False
    s = bisect.bisect_left(regions[chrom],start)
    e = bisect.bisect_right(regions[chrom],end)
    if s==e and s%2==0:
        return False
    else:
        if info_column:
            return tag[info_column]
        else:
            return True
    


def get_overlap_info_BEDtag2regions_df(tag,df,info_column=False):
    # if the file is too big, compare tags one by one
    # each regions[chr] is regionslist from regionslist_from_regions(regions)
    # return tag[column] if necessary
    chrom = tag[0]
    start = int(tag[1])
    end = int(tag[2])
    df = df[df['chrom']==chrom]
    # if df is a Empty DataFrame, the 'for' loop does not execute
    overlapped_id=[]
    for df_region in df.index:
        s = bisect.bisect_left([df.loc[df_region,'start'],df.loc[df_region,'end']],start)
        e = bisect.bisect_right([df.loc[df_region,'start'],df.loc[df_region,'end']],end)
        if s==e and s%2==0:
            pass
        else:
            overlapped_id.append(df_region)
    return overlapped_id
    
    
def get_overlap_info_regions2regions_strID(islands,regions):
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

def get_overlap_info_regions2regions_numID(islands,regions):
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
            regionStart,regionEnd = [],[]
            leftEnd=0
            for region in regions[chrom]:
                assert leftEnd <= region[0]<=region[1]            
                regionStart.append(region[0])
                regionEnd.append(region[1])
                leftEnd=region[1]
                            
            # for each island, check if overlapped with regions[chrom]
            for island in islands[chrom]:              
                assert island[0]<=island[1]
                e = bisect.bisect_left(regionEnd,island[0])
                s = bisect.bisect_right(regionStart,island[1])
                if s > e :
                    if chrom not in overlapped:
                        overlapped[chrom] = []
                    overlapped[chrom].append(island)
                else:
                    if chrom not in nonoverlapped:
                        nonoverlapped[chrom] = []
                    nonoverlapped[chrom].append(island)
                #i+=1
    return overlapped,nonoverlapped


def read_regions_from_bed_sep_strand(infile,species,expand=0,write_out=None,end_col=2,region_sorted=True):
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
                sline = line.strip().split()
                chrom = sline[0]
                start = int(sline[1])
                end = int(sline[end_col])
                strand = sline[5]
                if plus.match(strand):
                #if 1:
                    if chrom in chroms:
                        if chrom not in regions1:
                            regions1[chrom]=[]
                        regions1[chrom].append([max(0,start-expand),end+expand,'\t'.join(sline[3:])])
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



def read_regions_from_bed_not_sep_strand(infile,species,expand=0,write_out=False,end_col=2,region_sorted=True):
    if species in species_chroms.keys():
        chroms = species_chroms[species]
    else:
        print('Species not recognized!')
        exit(1)
    regions1,regions2 = {},{}
    with open(infile,'r') as inf:
        lines = inf.readlines(100000)
        while lines:
            #print(line)
            for line in lines:# and (if not re.match("#",lines)):
                #print(line)
                sline = line.strip().split()
                chrom = sline[0]
                start = int(sline[1])
                end = int(sline[end_col])
                #strand = sline[5]
                #if plus.match(strand):
                if 1:
                    if chrom in chroms:
                        if chrom not in regions1:
                            regions1[chrom]=[]
                        regions1[chrom].append([max(0,start-expand),end+expand,'\t'.join(sline[3:])])
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
            lines = inf.readlines(100000)
            
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
                        
    return regions1 



def read_regions_from_bed_by_chrom_not_sep_strand(infile,match_chrom,species,expand=0,write_out=None,end_col=2,region_sorted=True):
    if species in species_chroms.keys():
        chroms = species_chroms[species]
    else:
        print('Species not recognized!')
        exit(1)
    regions = []
    with open(infile,'r') as inf:
        i=0
        time1 = time.time()
        lines = inf.readlines(1024*1024)
        while lines:
            #print(line)
            for line in lines:
            #if not re.match("#",line):
                #print(line)
                sline = line.strip().split()
                chrom = sline[0]
                
                
                i+=1
                if i==10000000:
                    time2 = time.time()
                    print(time2-time1,'s')
                    i=0;time1 = time2
                
                #strand = sline[5]
                #if plus.match(strand):
                #if 1:
                if chrom == match_chrom:
                    start = int(sline[1])
                    end = int(sline[end_col])
                    regions.append([max(0,start-expand),end+expand,'\t'.join(sline[3:])])
            #print(regions1,regions2);exit(0)
            lines = inf.readlines(1024*1024)
            
    if region_sorted:
        #for chrom in regions:
        regions = sorted(regions,key = itemgetter(0),reverse=False)
    '''    
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
    '''                    
    return regions

def read_summits_from_narrowpeak_bed(infile,species,expand,write_out=None,end_col=2,region_sorted=True,peak_thre=4):
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
                if chrom in chroms and fold_enrichment>=peak_thre:
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

def fileter_BED_by_regions_not_sep_strand(infile,regions,species,fragmentsize=150,info_column=False,tagregions_sort = True):

    if species in species_chroms.keys():
        chroms = species_chroms[species]
    else:
        print('Species not recognized!')
        exit(1)
    filter_regions = {}
    shift = fragmentsize*0.5
    for chrom in regions:
        regions[chrom] = union_islands(regions[chrom],fragmentsize)
        regions[chrom] = regionslist_from_regions(regions[chrom])
    with open(infile,'r') as inf:
        lines = inf.readlines(1024*1024)
        while lines:
            for line in lines:
                #if not re.match("#",line):
                sline = line.strip().split()
                # if tag overlapped with regions, record tag info
                # all regions have been expanded in the union step, 
                # so no need to worry about the shifting of reads here
                if get_overlap_info_BEDtag2regions(sline,regions,info_column=False):
                    chrom = sline[0] 
                    start = int(sline[1])
                    end = int(sline[2]) 
                    strand = sline[5]                    
                    # for each chrom, build {} to save all tags
                    if chrom not in filter_regions:
                        filter_regions[chrom] = {}
                    # all start positions of plus strand with same end(key)
                    if plus.match(strand):
                        if end+shift not in filter_regions[chrom]:
                            filter_regions[chrom][end+shift]=set() # drop duplicate tags
                        filter_regions[chrom][end+shift].add(start+shift)
                    # all end positions of minus strand with same end(key)
                    if minus.match(strand):
                        #if start-shift not in filter_regions[chrom]:
                        #    filter_regions[chrom][start-shift]=set()
                        #filter_regions[chrom][start-shift].add(end-shift)  
                        if end-shift not in filter_regions[chrom]:
                            filter_regions[chrom][end-shift]=set()
                        filter_regions[chrom][end-shift].add(start-shift)
            lines = inf.readlines(1024*1024)
    # return [start,end] format tag_regions
    tag_regions = {}        
    for chrom in filter_regions:
        tag_regions[chrom] = [[start,end] for end in filter_regions[chrom] for start in filter_regions[chrom][end]]                   
        if tagregions_sort:
            tag_regions[chrom] = sorted(tag_regions[chrom],key = itemgetter(0))
    del filter_regions
    return tag_regions


 
def fileter_BED_by_regions_not_sep_strand_keep_score(infile,regions,species,fragmentsize=150,info_column=False,tagregions_sort = True):

    if species in species_chroms.keys():
        chroms = species_chroms[species]
    else:
        print('Species not recognized!')
        exit(1)
    filter_regions = {}
    shift = fragmentsize*0.5
    for chrom in regions:
        regions[chrom] = union_islands(regions[chrom],fragmentsize)
        regions[chrom] = regionslist_from_regions(regions[chrom])
    with open(infile,'r') as inf:
        lines = inf.readlines(1024*1024)
        while lines:
            for line in lines:
                #if not re.match("#",line):
                sline = line.strip().split()
                # if tag overlapped with regions, record tag info
                overlap_info = get_overlap_info_BEDtag2regions(sline,regions,info_column=info_column)
                if overlap_info:
                    chrom = sline[0] 
                    start = int(sline[1])
                    end = int(sline[2]) 
                    strand = sline[5]                    
                    # for each chrom, build {} to save all tags
                    if chrom not in filter_regions:
                        filter_regions[chrom] = {}
                    # all start positions of plus strand with same end(key)
                    if plus.match(strand):
                        if end+shift not in filter_regions[chrom]:
                            filter_regions[chrom][end+shift]=set() # drop duplicate tags
                        filter_regions[chrom][end+shift].add((start+shift,overlap_info))
                    # all end positions of minus strand with same end(key)
                    if minus.match(strand):
                        #if start-shift not in filter_regions[chrom]:
                        #    filter_regions[chrom][start-shift]=set()
                        #filter_regions[chrom][start-shift].add(end-shift)  
                        if end-shift not in filter_regions[chrom]:
                            filter_regions[chrom][end-shift]=set()
                        filter_regions[chrom][end-shift].add((start-shift,overlap_info))
                    #print(filter_regions);exit(0)
            lines = inf.readlines(1024*1024)
    # return [start,end] format tag_regions
    tag_regions = {}        
    for chrom in filter_regions:
        tag_regions[chrom] = [[start[0],end,start[1]] for end in filter_regions[chrom] for start in filter_regions[chrom][end]]                   
        #print(tag_regions);exit(0)
        if tagregions_sort:
            tag_regions[chrom] = sorted(tag_regions[chrom],key = itemgetter(0))
    del filter_regions
    return tag_regions
 
def write_regions_tab_join(outfile,regions1,regions2):

    with open(outfile,'w') as outf:
        for chrom in regions1:
            for ele in regions1[chrom]:
                outf.write('{}\t{}\n'.format(chrom,'\t'.join([str(i) for i in ele])))
        for chrom in regions2:
            for ele in regions2[chrom]:
                outf.write('{}\t{}\n'.format(chrom,'\t'.join([str(i) for i in ele])))

def write_regions_tab_join_not_sep_strand(outfile,regions1):

    with open(outfile,'w') as outf:
        for chrom in regions1:
            for ele in regions1[chrom]:
                outf.write('{}\t{}\n'.format(chrom,'\t'.join([str(i) for i in ele])))

    

def main(infile1,infile2,species,outfile1,outfile2):

    compare_regions1,compare_regions2 = read_region_from_bed(infile1,species)
    basic_regions1,basic_regions2 = read_region_from_bed(infile2,species)
    for ele in basic_regions1:
        basic_regions1[ele] = union_islands(basic_regions1[ele])
    for ele in basic_regions2:
        basic_regions2[ele] = union_islands(basic_regions2[ele])

    overlapped1,nonoverlapped1 = get_overlap_info_regions2regions_strID(compare_regions1,basic_regions1)        
    overlapped2,nonoverlapped2 = get_overlap_info_regions2regions_strID(compare_regions2,basic_regions2)
    write_regions_bed3(overlapped1,overlapped2,outfile1) 
    write_regions_bed3(nonoverlapped1,nonoverlapped2,outfile2)         
        
        
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--infile1', action = 'store', type = str,dest = 'infile1', help = 'input file to be compared/separated', metavar = '<file>')
    parser.add_argument('-b', '--infile2', action = 'store', type = str,dest = 'infile2', help = 'input file to be compared as basic', metavar = '<file>')
    parser.add_argument('-p','--outfile1', action = 'store', type = str,dest = 'outfile1', help = 'outfile of a overlapped with b', metavar = '<file>')
    parser.add_argument('-q','--outfile2', action = 'store', type = str,dest = 'outfile2', help = 'outfile of a not overlapped with b', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<9:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile1,args.infile2,args.species,args.outfile1,args.outfile2)
