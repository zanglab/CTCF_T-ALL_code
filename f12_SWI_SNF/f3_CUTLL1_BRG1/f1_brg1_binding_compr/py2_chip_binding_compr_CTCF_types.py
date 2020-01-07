import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
import CTCF_TALL_modules_new

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.4)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]




def window_smooth(df,half_window):
    smooth_df = pd.DataFrame(index=df.index,columns=df.columns)
    for col in np.arange(len(df.columns)):
        window_left = max(col-half_window,0)
        window_right = min(col+half_window,len(df.columns)-1)
        smooth_df.iloc[:,col] = df.iloc[:,window_left:window_right+1].mean(axis=1)
    return smooth_df   
    print(df,smooth_df);exit()
    return df



def plot_composite_figs(df,aml_gained,aml_lost,constitutive,outdir,basename,window):
    
    half_window = 3
    x = df.columns
    y1 = df.loc[aml_gained.index.intersection(df.index)];y1 = window_smooth(y1,half_window)
    y2 = df.loc[aml_lost.index.intersection(df.index)];y2 = window_smooth(y2,half_window)
#     y3 = df.loc[selected_ctrl];y3 = window_smooth(y3,half_window)
    y4 = df.loc[constitutive.index.intersection(df.index)];y4 = window_smooth(y4,half_window)
    print(y1.shape,y2.shape,y4.shape)
    fig = plt.figure(figsize=(2.8,2.8))
    plt.plot(x,y1.mean().values,label='T-ALL gained',color='red')
    plt.plot(x,y2.mean().values,label='T-ALL lost',color='blue')
#     plt.plot(x,y3.mean().values,label='Selected control',color='k')
    plt.plot(x,y4.mean().values,label='Constitutive',color='darkgrey')
#     sns.despine(offset=0, trim=False)
    plt.axes().set_xticks([x[0],x[int(len(x)/2)],x[-1]])
    plt.axes().set_xticklabels(['-{}kb'.format(str(window)[0]),0,'{}kb'.format(str(window)[0])],rotation=0)
    plt.title(basename.split(' ')[1],fontsize=16)
    plt.ylabel('ChIP-seq RPKM',fontsize=16)
    plt.legend(fontsize=14,bbox_to_anchor=[1.8,1.05],loc="upper right",borderaxespad=0.5,labelspacing=.5,handletextpad=0.5,handlelength=1,frameon=False)
    plt.savefig(outdir+os.sep+basename+'_{}.pdf'.format(window),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()
    


        
        
def main():

    indir = 'f1_chip_binding_pattern'
    outdir = 'f2_chip_binding_compr_CTCF_types'
    os.makedirs(outdir,exist_ok=True)
    
    aml_gained,aml_lost = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    constitutive_df = CTCF_TALL_modules_new.return_constitutive_df()
    
    
    ## composite plot
    for window in [1000]:
        suffix='_binding_pattern_w{}.csv'.format(window)
        infiles = glob.glob(indir+os.sep+'*{}'.format(suffix))
        for infile in infiles:
            basename = ' '.join(os.path.basename(infile).split(suffix)[0].split('_')[:3])
            print(basename)
            with open(infile) as inf:
                df = pd.read_csv(inf,sep='\t',index_col=0,header=None)
            plot_composite_figs(df,aml_gained,aml_lost,constitutive_df,outdir,basename,window)






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
