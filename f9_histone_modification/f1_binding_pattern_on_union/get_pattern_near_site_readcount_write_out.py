'''
Created on 04/15/2018
@authors: Chongzhi Zang, Zhenjia Wang<zw5j@virginia.edu>

This file is used to generate pattern collected from input tagsfile of bed/bam format 
for each start position width-expanded-site in input sitefile
'''

import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
#import association_with_regions
from get_reads_positions import reads_positions
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"

import re,bisect
plus = re.compile('\+')
minus = re.compile('\-')

    
def collect_reads_on_window(site_position,window,tag_regions):

    start,end = window[0],window[1]       

    if site_position.chr in tag_regions:
        s = bisect.bisect_left(tag_regions[site_position.chr],start)
        e = bisect.bisect_right(tag_regions[site_position.chr],end)
        return e-s
    return 0

def generate_patterns(outfile,site_positions,tagsfile,species,format,fragmentsize,width,bins,expand,mid):
    
    tags_regions,total = reads_positions(tagsfile,species,format,fragmentsize)
    half_window_size = int(width/bins)
    #sum_pattern = pd.DataFrame(columns = np.arange(-1*width,width+1,2*width/bins,dtype=int))
    outf = open(outfile,'w')
    outf.write('id\t{}\n'.format('\t'.join(map(str,np.arange(-1*width,width+1,2*width/bins,dtype=int)))))
    for position_index in site_positions.index:
        #if position_index%1000==0:
        #    print(position_index)
        site_position = site_positions.loc[position_index] 
        pattern = [] 
        if mid:
            point_position = int((site_position.position + site_position.end)/2)
        else:
            point_position = site_position.position
        #print(mid,point_position,site_position.position,site_position.end);exit(0)  
        for ii in np.arange(point_position-width,point_position+width+1,2*width/bins,dtype=int):
            window = [ii-half_window_size-expand,ii+half_window_size+expand]#;print(window)
            tag_nums = collect_reads_on_window(site_position,window,tags_regions)
            pattern.append(round(tag_nums*1000000000/(total*2*(half_window_size+expand)),2))
        if minus.match(site_position.strand):
            pattern.reverse()
        outf.write('{}\t{}\n'.format(site_position.id,'\t'.join(map(str,pattern))))
        #sum_pattern.loc[site_position.id] = pattern
        #plt.figure()
        #plt.plot(pattern)
        #plt.savefig('{}_{}.png'.format(site_position.id,site_position.strand))
        #plt.close()
    #sum_pattern = sum_pattern.round(2)
    #return sum_pattern
    outf.close()
    return 1

def main(sitefile,tagsfile,outfile,species,format,fragmentsize,width,bins,expand,mid):


    # site_positions.columns should at least include: ['chr','position','strand','id']
    site_positions = pd.read_csv(args.sitefile,header = None,sep='\t')
    site_positions.columns = ['chr','position','end','id','val','strand']+[i for i in np.arange(site_positions.shape[1]-6)]
    patterns = generate_patterns(outfile,site_positions,tagsfile,species,format,fragmentsize,width,bins,expand,mid)
    #patterns.to_csv(outfile,sep='\t')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--sitefile', action = 'store', type = str,dest = 'sitefile', help = 'input bed file of sites', metavar = '<file>')
    parser.add_argument('-t', '--tagsfile', action = 'store', type = str,dest = 'tagsfile', help = 'input bed/bam file of reads/tags', metavar = '<file>')
    parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of reads pattern for each site', metavar = '<file>',required=True)
    parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    parser.add_argument('-f', '--format', action = 'store', type = str,dest = 'format', help = 'format of file, BED or BAM', metavar = '<str>',required=True)
    parser.add_argument('-g', '--fragmentsize', action = 'store', type = int,dest = 'fragmentsize', help = 'fragmentsize for the shift of reads. Default: 147', metavar = '<int>',default=147)    
    parser.add_argument('-w','--width', action = 'store', type = int,dest = 'width', help = 'one side width for the pattern plot, default: 1000', metavar = '<int>',default=1000)
    parser.add_argument('-b','--bins', action = 'store', type = int,dest = 'bins', help = 'bins to separate the width, default: 100', metavar = '<int>',default=100)
    parser.add_argument('-e', '--expand', action = 'store', type = int,dest = 'expand', help = 'one side expand region of each bin. Default: 0', metavar = '<int>',default=0)
    parser.add_argument('-m', '--mid', action = 'store_true', dest = 'mid', help = 'whether to use middle side for expansion. Default: False',default=False)
    
    args = parser.parse_args()
    if(len(sys.argv))<11:
        parser.print_help()
        sys.exit(1)
  
    main(args.sitefile,args.tagsfile,args.outfile,args.species,args.format,args.fragmentsize,args.width,args.bins,args.expand,args.mid)
