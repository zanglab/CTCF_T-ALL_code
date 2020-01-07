import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
sns.despine(offset=0, trim=True)
import CTCF_TALL_modules


def get_mutation_info(infile):
    
    with open(infile) as inf:
        df = pd.read_csv(inf,sep='\t',index_col=0,header=None)
    avg_mutation_rate = [i/df.shape[0] for i in df.sum()]    
    return avg_mutation_rate
        

def composite_plot(plot_lists,plot_colors,figname):

    plt.figure()
    for ii in np.arange(len(plot_lists)):
        plot_list = plot_lists[ii]
        x = np.arange(len(plot_list))#;print(plot_list,x)
        color = plot_colors[ii]
        plt.plot(x,plot_list,color)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600)
    plt.close()

def main():

    outdir = 'f3_mutation_compr'
    os.makedirs(outdir,exist_ok=True)
    
    for celltype in ['CUTLL1','Jurkat','PD30','PD9']:
        for expand in [9,200]:
            percentage_df = pd.DataFrame(columns = np.arange(1,2*expand+2))
            gained_file ='f2_mutation_on_cancer_specific_CTCF/{}_mutation_on_TALL_{}_expand{}.mutation_rate.csv'.format(celltype,'gained',expand)
            lost_file =  'f2_mutation_on_cancer_specific_CTCF/{}_mutation_on_TALL_{}_expand{}.mutation_rate.csv'.format(celltype,'lost',expand)
            union_file = 'f1_mutation_around_CTCF/{}_mutation_on_CTCF_expand{}.csv.mutation_event_rate.csv'.format(celltype,expand)
            
            gained_mutation_avg_rate = get_mutation_info(gained_file)
            lost_mutation_avg_rate = get_mutation_info(lost_file)
            union_mutation_avg_rate = get_mutation_info(union_file)
            
            percentage_df.loc['union'] = union_mutation_avg_rate
            percentage_df.loc['lost'] = lost_mutation_avg_rate
            percentage_df.loc['gained'] = gained_mutation_avg_rate
            percentage_df.to_csv(outdir+os.sep+'{}_expand{}.csv'.format(celltype,expand))
            
#             plot_lists=[gained_percentage,lost_percentage,const_percentage,union_percentage]
#             plot_colors = ['r','b','k','grey']
#             plot_lists=[gained_percentage,lost_percentage,union_mutation_percentage]
#             plot_colors = ['r','b','grey',]
#             figname = outdir+os.sep+'{}_expand{}.png'.format(celltype,expand)
            #composite_plot(plot_lists,plot_colors,figname)
            #exit()
            #print(const_percentage)
            
        
        
        

            
            
            
            


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
