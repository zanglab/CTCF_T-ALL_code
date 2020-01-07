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
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")

#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"

def scatter_plot(df,marked_index,figname,mark_col):
    
    pylist = df.index
    plt.figure(figsize=(2.8,2.8))
    for index in df.index:
        xp = list(df.index).index(index)
        values = -1*np.log10(df.loc[index,mark_col])#;print(values)
        plt.scatter(xp,values,color='k',s=6)
    
    max_p = -1*np.log10(df.iloc[0,-1])
    values_reset = max_p*1.05
    # top 3 and TFs in shared top 10
    sorted_marker_ids = np.append(df.index[:2],df.index[:20].intersection(marked_index)[:5])
    sorted_marker_pos = sorted(set([list(df.index).index(marker) for marker in  sorted_marker_ids]))
    x_pos_m=0
    for marker_id in sorted_marker_pos:
        xp = marker_id #;print(xp,df)
        values = -1*np.log10(df.loc[df.index[marker_id],mark_col])#;print(values)
        plt.scatter(xp,values,color='red',s=15)
        # ==== mark the index label, not overlap with each other
        values_reset = min(values,values_reset-max_p*0.1)#;print(values,values_reset)
        # ==== mark the index and plot the arrow separately
        plt.text(xp+120+x_pos_m,values_reset,df.index[marker_id])
        if df.index[marker_id] in marked_index:
            plt.text(xp+120+x_pos_m,values_reset,df.index[marker_id],color='red')
        # ==== plt.arrow(x,y,dx,dy)
        plt.arrow(xp+120+x_pos_m,values_reset+max_p*0.03,-90-x_pos_m,values-values_reset-max_p*0.03,\
                  length_includes_head = True,head_width=max_p*0.02,head_length=35,fc='k',ec='k')
        x_pos_m = x_pos_m + 30

    #plt.title('Predicted co-factors of \nT-ALL gained bindings')
    plt.xlabel('TF Rank',fontsize=16)
    plt.ylabel('-log$_{{10}}$ $P$ value',fontsize=16)
    plt.axes().set_xticks([0,len(df.index)])
    plt.axes().set_xticklabels([1,len(df.index)],rotation=0, ha='center',fontsize=16,color='k')
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()




def main(infile):

    outdir = 'rank_dot_figs'
    os.makedirs(outdir,exist_ok=True)
    
    marked_index = ['NOTCH1']
#     infile = '../bart_output/PAML_14_combined_DNAme_gene_removed_bart_results.txt'
    basename = os.path.basename(infile).split('.')[0]
    df = pd.read_csv(infile,index_col=0,sep=',')
    mark_col = df.columns[-1]
    df = df.sort_values(by=[mark_col],ascending=True)
    figname = outdir+os.sep+basename+'.pdf'
    scatter_plot(df,marked_index,figname,mark_col)

    
    





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
