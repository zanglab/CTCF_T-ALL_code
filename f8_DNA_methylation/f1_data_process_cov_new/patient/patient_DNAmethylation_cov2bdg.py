import sys,argparse
import os,glob
import numpy as np
import pandas as pd
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

def plot_methy_count(df,basename):
    '''plot the distribution of all count for all sites'''
    total = df['methy_count']+df['methy_count']
    total = [i for i in total.values if i<=100]
    plt.figure()
    plt.hist(total,bins=100)
    plt.title(basename)
    plt.xlabel('count')
    plt.yscale('log')
    plt.xlim([0,100])
    plt.savefig('methy_count_fig/{}.png'.format(basename),bbox_inches='tight',pad_inches=0.1,transparent=True)


def cov2bdg(df,basename,count_filter):

    df['total'] = df['methy_count']+df['unmethy_count']
    df = df[df['total']>=count_filter]
    df['chr'] = 'chr'+df['chr']
    df = df[['chr','start','end','methy_percentage']]
    df.to_csv('patient_bdg/{}_countfilter_{}.bdg'.format(basename,count_filter),header=None,index=False,sep='\t')
    #print(df);exit()

def main():

    methylationfiles = glob.glob('/nv/vol190/zanglab/zw5j/projects_data/T_ALL_CTCF_panos/panos/DNA-methylation/Patient/*cov')
    for methylationfile in methylationfiles:
        basename = os.path.basename(methylationfile).split('_')[1]
        with open(methylationfile) as inf:
            df = pd.read_csv(inf,sep='\t',header=None,low_memory=False)
        df.columns = ['chr','start','end','methy_percentage','methy_count','unmethy_count']
        #plot_methy_count(df,basename)
        
        for count_filter in [5,10]:
            cov2bdg(df,basename,count_filter)
        
        




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
  
    main()
