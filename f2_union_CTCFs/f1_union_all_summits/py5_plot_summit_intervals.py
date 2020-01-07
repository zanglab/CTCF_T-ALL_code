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




def plt_pdf(mylist,bins=0,ylog=False,xlog=False,figname=None,bar=False):
    mylist = sorted(mylist)
    a = np.min(mylist)
    b = np.max(mylist)
    if bins!=0:
        sep=(b-a)/bins
    else:
        sep = (b-a)/len(mylist) 
    x = np.arange(a,b+sep,sep)
    y=[]
    l=0
    for xi in x:
        r = bisect.bisect_right(mylist,xi)
        y.append(r-l)
        l=r
    #y.append(len(mylist)-l)
    #plt.figure()
    if bar:
        plt.bar(x,y,width=x[1]-x[0])
    else:
        if ylog:
            plt.yscale('log')
        if xlog:
            plt.xscale('log')
        #plt.xlim([0,500])
        #plt.plot(x,y)
    if figname:
        plt.savefig(figname,bbox_inches='tight',pad_inches=0.02)
    #plt.close()
    return x,y



def main(infile):               

    subintervals = [int(line.strip()) for line in open('f1_peak2000_datasets_union_summits/all_summits_fe4_intervals.txt').readlines()]
    subintervals = [i for i in subintervals if i<500]

   
    plt.figure(figsize=(3,3))
    x,y=plt_pdf(subintervals,bins=100)
    plt.plot(x,y)
    #g = sns.distplot(subintervals,kde=False)  
    plt.yscale('log')
    # plt.ylabel('#')
    plt.savefig('summits_intervals.pdf',bbox_inches='tight',pad_inches=0.02)
    plt.close()





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