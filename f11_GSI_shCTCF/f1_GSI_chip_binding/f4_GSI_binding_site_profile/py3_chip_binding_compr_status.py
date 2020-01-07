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




def window_smooth(df,half_window):
    smooth_df = pd.DataFrame(index=df.index,columns=df.columns)
    for col in np.arange(len(df.columns)):
        window_left = max(col-half_window,0)
        window_right = min(col+half_window,len(df.columns)-1)
        smooth_df.iloc[:,col] = df.iloc[:,window_left:window_right+1].mean(axis=1)
    return smooth_df   
    print(df,smooth_df);exit()
    return df



def plot_composite_figs(df1,df2,aml_gained,constitutive,outdir,basename,window):
    
    half_window = 3
    x = df1.columns
    y11 = df1.loc[aml_gained.index.intersection(df1.index)];y11 = window_smooth(y11,half_window)
    y21 = df1.loc[constitutive.index.intersection(df1.index)];y21 = window_smooth(y21,half_window)

    y12 = df2.loc[aml_gained.index.intersection(df2.index)];y12 = window_smooth(y12,half_window)
    y22 = df2.loc[constitutive.index.intersection(df2.index)];y22 = window_smooth(y22,half_window)

    # print(y1.shape,y2.shape,y4.shape)
    fig = plt.figure(figsize=(4,3))
    plt.plot(x,y11.mean().values,label='T-ALL gained {}'.format(basename.split(' ')[0]),color='red')
    plt.plot(x,y21.mean().values,label='Constitutive {}'.format(basename.split(' ')[0]),color='darkgrey')

    plt.plot(x,y12.mean().values,label='T-ALL gained {}'.format(basename.split(' ')[-1]),color='red',ls='--')
    plt.plot(x,y22.mean().values,label='Constitutive {}'.format(basename.split(' ')[-1]),color='darkgrey',ls='--')

    sns.despine(offset=0, trim=False)
    plt.axes().set_xticks([x[0],x[-1]])
    plt.axes().set_xticklabels(['-{}kb'.format(str(window)[0]),'{}kb'.format(str(window)[0])],rotation=30)
    plt.title(basename,fontsize=18)
    plt.ylabel('RPKM',fontsize=18)
    plt.legend(fontsize=14,bbox_to_anchor=[1.5,1.0],loc="upper right",borderaxespad=0.5,labelspacing=.5,handletextpad=0.5,handlelength=1,frameon=False)
    plt.savefig(outdir+os.sep+basename+'_{}.pdf'.format(window),bbox_inches='tight',pad_inches=0.1,dpi=600)
    plt.close()
    


        
        
def main():

    indir = 'f1_chip_binding_pattern'
    outdir = 'f3_chip_binding_compr_status'
    os.makedirs(outdir,exist_ok=True)
    
    aml_gained,aml_lost = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    constitutive_df = CTCF_TALL_modules_new.return_constitutive_df()
    

    ## composite plot comparing DMSO vs. GSI
    infile1='f1_chip_binding_pattern/dmso_binding_pattern_w1000.csv'
    infile2='f1_chip_binding_pattern/gsi_3d_binding_pattern_w1000.csv'
    basename='DMSO vs. GSI'
    for window in [1000]:
        df1 = pd.read_csv(infile1,sep='\t',index_col=0,header=None)
        df2 = pd.read_csv(infile2,sep='\t',index_col=0,header=None)
        plot_composite_figs(df1,df2,aml_gained,constitutive_df,outdir,basename,window)


    infile1='f1_chip_binding_pattern/gsi_3d_binding_pattern_w1000.csv'
    infile2='f1_chip_binding_pattern/gsi_3d_w4h_binding_pattern_w1000.csv'
    basename='GSI vs. GSI-wo'
    for window in [1000]:
        df1 = pd.read_csv(infile1,sep='\t',index_col=0,header=None)
        df2 = pd.read_csv(infile2,sep='\t',index_col=0,header=None)
        plot_composite_figs(df1,df2,aml_gained,constitutive_df,outdir,basename,window)


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
