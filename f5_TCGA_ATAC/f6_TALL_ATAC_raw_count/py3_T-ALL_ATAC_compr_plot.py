import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
#import association_with_regions
from get_reads_positions import reads_positions
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
from scipy import stats
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False,'grid.color': 'grey'})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
import CTCF_TALL_modules_new
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset']='custom'
matplotlib.rcParams['mathtext.rm']='Arial'


def stats_test_on_signals(df,cancer_columns,ctrl_columns):
    cancer_df = df[cancer_columns]
    ctrl_df = df[ctrl_columns]
    s,p = stats.ttest_ind(cancer_df,ctrl_df,axis=1)
    s = np.nan_to_num(s)#;print(s,p);exit()
    fc = cancer_df.mean(axis=1)-ctrl_df.mean(axis=1)
    return fc,s,p


def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),99.7)*1.3 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),.2)*1.3
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label = '{:.1e}'.format(p)
    p_label='*'
    if p<0.001:
        p_label = '**'

    if p<0.05:
        if p_label[-2]=='0':
            p_label = p_label[:-2]+p_label[-1]
        if compr_pos[2] == 't':
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*h, p_label, ha='center', va='bottom', color=col,fontsize=15) # "{:.1e}".format(p)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*1.1, y2*1.1, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*1.8, p_label, ha='center', va='bottom', color=col,fontsize=15)

   
def main():

    signal_df = pd.read_csv('f2_combined_raw_count_csv/Jurkat_CD4_ATAC_combined_raw_count_log2_QN.csv',index_col=0)

    
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    const_df = CTCF_TALL_modules_new.return_constitutive_df()
    
    const_df = signal_df.loc[const_df.index]
    gained_df = signal_df.loc[gained_df.index]
    lost_df = signal_df.loc[lost_df.index]
    #print(gained_df);exit()
    
    ctrl_columns = ['GSM1155964','GSM1155965','GSM1155966','GSM1155967','GSM1155968','GSM1155969'] # CD4
    cancer_columns = ['GSM2411156','GSM2411157','GSM2411158'] # Jurkat
     
    figname = 'Jurkat_vs_CD4_ATAC_signal_on_TALL_CTCF.pdf'
    plt.figure(figsize=(2.4,2.4))
    gained_fc,gained_s,gained_p = stats_test_on_signals(gained_df,cancer_columns,ctrl_columns)
    lost_fc,lost_s,lost_p = stats_test_on_signals(lost_df,cancer_columns,ctrl_columns)
    const_fc,const_s,const_p = stats_test_on_signals(const_df,cancer_columns,ctrl_columns)
#     box_vals = [lost_s,const_s,gained_s]
    box_vals = [lost_fc,const_fc,gained_fc]
    positions = [1,2,3]
    g = plt.boxplot(box_vals,positions=positions,widths = .55,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='k'),showfliers=False)
    

    for compr_pos in [[0,1,'b'],[1,2,'t']]:
        mark_pvalue(compr_pos,positions,box_vals)
    colors = [ 'blue','darkgrey','red']
#     for patch, color in zip(g['boxes'], colors):
#         patch.set_facecolor(color)
#         patch.set_alpha(0.9)
    scatter_X = []
    for position_id in np.arange(len(positions)):
        scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
        plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=20,zorder=0,alpha=1,marker='o')
    
    
    plt.axhline(y=0,color='grey',linestyle='--',linewidth=.7)
#    plt.ylabel('ATAC-seq \ndifferential score',fontsize=18)
    plt.ylabel('$\Delta$ (ATAC-seq score)\n T-ALL over T-cell')
#     plt.axes().tick_params(axis='y',direction='out', length=4, width=.8, colors='black')
#     plt.axes().tick_params(axis='x',direction='out', length=4, width=.8, colors='black')
#     plt.ylim([-5,4])

#     plt.title('T-ALL',fontsize=18)
#     sns.despine(offset=None, trim=False)
    plt.axes().set_xticklabels(['T-ALL lost','Constitutive','T-ALL gained'],rotation=30,ha='center',fontsize=16) 
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
