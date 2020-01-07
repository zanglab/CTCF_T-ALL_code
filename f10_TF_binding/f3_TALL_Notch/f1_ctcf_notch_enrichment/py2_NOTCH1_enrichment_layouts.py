import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
from scipy import stats
import CTCF_TALL_modules_new
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]

def mark_p(p,position,val):
    if p<0.05:
        star_mark="*"
    if p<0.001:
        star_mark="**"
    if p<0.05:
        plt.text(position,val-2,star_mark,ha='center',fontsize=17)
#     if 0.0001<=p<0.001:
#         plt.text(position,val,'**',ha='center')
#     if p<0.0001:
#         plt.text(position,val,'***',ha='center')



def bar_plot_enrichment(gs,ls,cs,figname):
#     a0,a1,a2,a3,a4,a5 = val_gained
#     b0,b1,b2,b3,b4,b5 = val_lost
#     c0,c1,c2,c3,c4,c5 = val_const
    xticklabels = ['NOTCH1','NOTCH1-D','MYC','NOTCH1 & MYC','NOTCH1-D & MYC']
    xticklabels = ['w/ NOTCH1','w/ dyNOTCH1','w/ MYC']
    xticklabels = ['Containing NOTCH1','Containing\ndynamic NOTCH1']
    plt.figure(figsize=(3,3))
    for ii in np.arange(1,len(gs)-3): 
        plot_positions = [ii+0.1,ii+0.35,ii+0.6]
        plot_vals = [100*cs[ii]/cs[0],100*ls[ii]/ls[0],100*gs[ii]/gs[0]]
        g=plt.bar(plot_positions,plot_vals,width=0.22,color = ['silver','blue','red'],linewidth=0)
        s,p = stats.fisher_exact([[cs[0]-cs[ii],cs[ii]],[ls[0]-ls[ii],ls[ii]]]);print(xticklabels[ii-1],'lost:',s,p)
        mark_p(p,ii+0.35,100*ls[ii]/ls[0])
        s,p = stats.fisher_exact([[cs[0]-cs[ii],cs[ii]],[gs[0]-gs[ii],gs[ii]]]);print(xticklabels[ii-1],'gained:',s,p)
        mark_p(p,ii+0.6,100*gs[ii]/gs[0])
    plt.axes().set_xticks(np.arange(1,len(gs)-2)+0.5)
    plt.ylabel('CTCF-containing\n domains (%)',fontsize=16)
    legends = ['Constitutive','T-ALL lost','T-ALL gained']
    plt.legend((g[0],g[1],g[2]),legends,fontsize=14,loc='upper left', bbox_to_anchor=[0.0,1.2,.2,.1],frameon=False,borderaxespad=0.0,labelspacing=.1,handletextpad=0.1,handlelength=1)
    plt.xlim([.8,2.9])
    plt.ylim([0,85])
    #plt.axes().tick_params(axis='x',direction='out', length=3, width=.8, colors='black')
#     sns.despine(offset=0, trim=False)
    plt.axes().tick_params(axis='x',direction='out', length=0, width=.8, colors='black')
    plt.axes().set_xticklabels(xticklabels,rotation=20,ha = 'right',fontsize=16)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.02,dpi=600,transparent=True)
    plt.close()



def return_co_localization_info(binding_df,flag):

    df = binding_df[['if_intra_domain_notch', 'if_intra_domain_dynamic_notch','if_intra_domain_MYC']]
    df.columns = ['N','DN','M']    
    a0 = df.shape[0]
    a1 = df[(df['N']==1)].shape[0]  
    a2 = df[(df['DN']==1)].shape[0] 
    a3 = df[(df['M']==1)].shape[0] 
    a4 = df[(df['N']==1)&(df['M']==1)].shape[0] 
    a5 = df[(df['DN']==1)&(df['M']==1)].shape[0] 
    
    print(flag,'\tTotal:\t',a0)
    print(flag,'\tNotch:\t',a1)
    print(flag,'\tNotch-d:\t:',a2)
    print(flag,'\tMYC:\t',a3)
    print(flag,'\tNotch-Myc:\t',a4)
    print(flag,'\tNotch-d-Myc:\t',a5);print()
    return a0,a1,a2,a3,a4,a5


def main(infile):

    outdir = 'f2_figs_notch_enrichment_layouts'
    os.makedirs(outdir,exist_ok=True) 

    union_feature_df = CTCF_TALL_modules_new.return_cancer_specific_combined_features('T-ALL')
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    constitutive_df = CTCF_TALL_modules_new.return_constitutive_df()
    
    # cols = 'if_intra_domain_notch', 'if_intra_domain_dynamic_notch','if_intra_domain_MYC'
    val_gained = return_co_localization_info(union_feature_df.loc[gained_df.index],'gained')
    val_lost = return_co_localization_info(union_feature_df.loc[lost_df.index],'lost')
    val_const = return_co_localization_info(union_feature_df.loc[constitutive_df.index],'const')
    
    figname = '{}/NOTCH1_enrichment_bar.pdf'.format(outdir)
    bar_plot_enrichment(val_gained,val_lost,val_const,figname)



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
