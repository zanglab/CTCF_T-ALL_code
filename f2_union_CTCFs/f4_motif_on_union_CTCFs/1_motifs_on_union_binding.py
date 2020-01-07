import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import association_with_regions
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})

import re,bisect
plus = re.compile('\+')
minus = re.compile('\-')
import CTCF_TALL_modules_new



def main(infile):

    '''this is used to find the motifs on each of the union bindings,
    for each binding site, only  the chr,start,end,pvalue,strand position of the motif with smallest
    p value was kept, and the 6th column is the number of motifs within this region'''
   # if species in species_chroms.keys():
   #     chroms = species_chroms[species]
    
    species = 'hg38'

    outdir='f1_union_binding_motif'
    os.makedirs(outdir,exist_ok=True)
    
#     collection_df = CTCF_TALL_modules_new.return_collection_df()
    union_df = CTCF_TALL_modules_new.return_union_binding_bed()
    
    motif_file = '/nv/vol190/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/CTCF_MA0139.1.bed'
    motif_regions_plus,motif_regions_minus = association_with_regions.read_regions_from_bed_sep_strand(motif_file,species)
    startlist_plus,startlist_minus = {},{}
    for ele in motif_regions_plus:
        startlist_plus[ele] = [i[0] for i in motif_regions_plus[ele]];#print(motif_regions_plus);exit(0)
    for ele in motif_regions_minus:
        startlist_minus[ele] = [i[0] for i in motif_regions_minus[ele]]

    bedfile = open('{}/union_CTCF_with_motif.bed'.format(outdir),'w')
    bedfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('chr','start','end','id','motif_start','motif_end','motif_pvalue','motif_strand','#motifs'))
    for binding_id in union_df.index:
        overlapped = []
        binding = union_df.loc[binding_id]
        chr = binding.chr
        start = int(binding.start)
        end = int(binding.end)
        if chr in startlist_plus:
            s = bisect.bisect_left(startlist_plus[chr],start-19)
            e = bisect.bisect_right(startlist_plus[chr],end)
            if s!=e:
                for i in np.arange(s,e):
                    overlapped.append(motif_regions_plus[chr][i])

        # minus strand
        if chr in startlist_minus:
            s = bisect.bisect_left(startlist_minus[chr],start-19) # end list vs. start == start list vs. start-19
            e = bisect.bisect_right(startlist_minus[chr],end)
            if s!=e:
                for i in np.arange(s,e):
                    overlapped.append(motif_regions_minus[chr][i])
        '''find the motif with smallest pvalue'''
        if len(overlapped)!=0:
            overlapped = sorted(overlapped,key = lambda x: float(x[2].split('\t')[1]),reverse=False)
            motif = overlapped[0]
            bedfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr,start,end,binding_id,motif[0],motif[1],motif[2].split('\t')[1],motif[2].split('\t')[2],len(overlapped)))
        else:
            bedfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr,start,end,binding_id,'N','N','N','N','N'))
    
#         if binding_id%10000==0:
#             print(binding_id)
            
    bedfile.close()




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
  
    main(args.infile)
