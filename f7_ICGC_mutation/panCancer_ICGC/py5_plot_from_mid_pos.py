import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
# matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
# matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
sns.despine(offset=0, trim=True)
# import CTCF_TALL_modules
from scipy import stats
matplotlib.rcParams["font.sans-serif"] = ["Arial"]

        
def win_smooth_plot(df,index,color,pos_add,window=40,smooth_window=0):
    mid_pos = int(len(df.columns)/2)
    x_range  = np.arange(mid_pos,len(df.columns)-1,window)
    y_list = []
    x_ticklabels=[]
    for x in x_range:
        mirror_x = len(df.columns)-x-1
        x_ticklabels.append('[{},{})'.format(x-smooth_window-mid_pos,min(len(df.columns)-1-mid_pos,x+window+smooth_window-mid_pos)))
        plot_window1 = np.arange(x-smooth_window,min(len(df.columns)-1,x+window+smooth_window))
        plot_window2 = np.arange(max(0,mirror_x-window-smooth_window),mirror_x+smooth_window)
        plot_col = df.columns[np.append(plot_window1,plot_window2)]
        y_list.append(df.loc[index,plot_col].mean())
#    g = plt.plot(x_range,y_list,color=color)
    g = plt.bar(x_range+window*pos_add,y_list,color=color,width=window*.25,lw=0)
    x_ticks = x_range+window*pos_add
    return g,x_ticks,x_ticklabels     
    
        

def composite_plot(df,figname,window,title):
    plt.figure(figsize=(6,5))
    g1,x_ticks,x_ticklabels = win_smooth_plot(df,'union','grey',0,window)
    g2,x_ticks2,x_ticklabels = win_smooth_plot(df,'lost','b',.25,window)
    g3,x_ticks3,x_ticklabels = win_smooth_plot(df,'gained','r',.5,window)
    plt.legend([g1[0],g2[0],g3[0]],['Union','Lost','Gained'],loc='upper right',frameon=False,fontsize=18,borderaxespad=-0.1,labelspacing=.05,handlelength=1,handletextpad=0.1)
    mid = int(len(df.columns)/2)
    plt.axes().set_xticks(x_ticks3)
    plt.axes().set_xticklabels(x_ticklabels,rotation=30,ha='right')
    plt.ylabel('mutation rate')
    plt.xlabel('distance to motif center (bp)')
    plt.ylim([0,0.000006])
    plt.title(title)
    plt.axes().tick_params(axis='x',direction='out', length=0, width=1, colors='black')
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600)
    plt.close()


def main():

    indir = 'f3_mutation_compr'
    outdir = 'f5_from_motif_mid_pos_replot_f3_figs'
    os.makedirs(outdir,exist_ok=True)
    cancertype_mutation_matchness = {'BRCA':'BRCA','CRC':'COAD','LUAD':'LUAD','PRAD':'PRAD','AML':'AML','PRAD_TissueAdded':'PRAD'}
    cancertypes=['BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']
    
    for celltype in cancertypes:
        for expand in [200]:
            percentage_df = pd.read_csv(indir+os.sep+'{}_expand{}.csv'.format(celltype,expand),index_col=0)
            
            # paired t-test on mutation rate
            print('\n{}, expand {}'.format(celltype,expand))
            s,p = stats.ttest_rel(percentage_df.loc['gained'],percentage_df.loc['union'])
            print('gained, s\t{}'.format(s))
            print('gained, p\t{}'.format(p))
            s,p = stats.ttest_rel(percentage_df.loc['lost'],percentage_df.loc['union'])
            print('lost, s\t{}'.format(s))
            print('lost, p\t{}'.format(p))

               
            for window in [40]:
                if window<expand and (window > expand/6.0):
                    figname = outdir+os.sep+'{}_expand{}_window{}.pdf'.format(celltype,expand,window)
                    composite_plot(percentage_df,figname,window,cancertype_mutation_matchness[celltype])
                    
        
        
        

            
            
            
            


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
