'''
This file is used to compr the CTCF binding signal 
in cancer cell lines
between gained/lost/ctrl for each cancer type
'''

import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
#import association_with_regions
#from get_reads_positions import reads_positions
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
from scipy import stats
import seaborn as sns
sns.set(font_scale=1.6)
sns.set_style("whitegrid", {'axes.grid' : False,'grid.color': 'grey'})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
import CTCF_TALL_modules



def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit');print(s,p)
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),99.5)*1.1,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),1)*1.1
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='{:.1e}'.format(p)
#     if p_label[-2]=='0':
#         p_label = p_label[:-2]+p_label[-1]
    if p<0.05:

        if compr_pos[2] == 't':
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*h,  p_label, ha='center', va='bottom', color=col,fontsize=16) # "{:.1e}".format(p)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*1.1, y2*1.1, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*2, p_label, ha='center', va='bottom', color=col,fontsize=16)

def signal_compr(box_vals,figname,title):
        
    # plot to compare binding signals
    plt.figure(figsize=(2,2.7))
    positions = np.arange(len(box_vals))
    colors = ['silver','red']
    g = plt.boxplot(box_vals,positions=positions,widths = .55,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='k'),showfliers=False)
    for compr_pos in [[0,1,'t',1.1]]:
        mark_pvalue(compr_pos,positions,box_vals)

    # scatters    
    scatter_X = []
    for position_id in np.arange(len(positions)):
        scatter_x = np.random.normal(positions[position_id],0.06,len(box_vals[position_id]))
        plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=20,zorder=0,alpha=1,marker='o')
    
    #plt.axhline(y=0,color='grey',linestyle='--',linewidth=.7)
    plt.ylabel('CTCF binding levels',fontsize=18)
#     plt.axes().tick_params(axis='y',direction='out', length=4, width=.8, colors='black')
#     plt.axes().tick_params(axis='x',direction='out', length=4, width=.8, colors='black')
    plt.ylim([-3,80])
    if title =="PRAD_TissueAdded":
        title="PRAD"
    plt.title(title,fontsize=20)
    sns.despine(offset=None, trim=False)
    plt.axes().set_xticklabels(['Constitutive','Gained'],rotation=30,ha='right',fontsize=18) 
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()


   
def main():

    indir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f5_binding_signal_compr/f3_gained_constitutive_CTCF_binding_signal_compr'
    outdir='figures/cancer_specific_binding_signal_compr_panCancer_gained_constitutive'
    os.makedirs(outdir,exist_ok=True)

#     cancertype_title = {'T-ALL':'T-ALL', 'Lung_cancer':'LUAD', 'Colon_cancer':'CRC', 'Breast_cancer':'BRCA'}
    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']
    for cancerType in cancertypes:
        gained_signals = np.load(indir+os.sep+'binding_level_on_cancer_datasets_{}.gained_signals.npy'.format(cancerType))
        lost_signals = np.load(indir+os.sep+'binding_level_on_cancer_datasets_{}.lost_signals.npy'.format(cancerType))
        constitutive_signals = np.load(indir+os.sep+'binding_level_on_cancer_datasets_{}.constitutive_signals.npy'.format(cancerType))
        
#         box_vals = [lost_signals,constitutive_signals,gained_signals]
        box_vals = [constitutive_signals,gained_signals]
#         title = cancertype_title[cancerType]
        title = cancerType
        figname = outdir+os.sep+'gained_constitutive_binding_compr_{}.png'.format(cancerType)
        signal_compr(box_vals,figname,title)






if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()

