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
import CTCF_TALL_modules_new

def stats_test_on_signals(df,cancer_columns,ctrl_columns):
    cancer_df = df[cancer_columns]
    ctrl_df = df[ctrl_columns]
    s,p = stats.ttest_ind(cancer_df,ctrl_df,axis=1)
    s = np.nan_to_num(s)#;print(s,p);exit()
    fc = cancer_df.mean(axis=1)-ctrl_df.mean(axis=1)
    return fc,s,p


def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),100)*compr_pos[3] ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),1)*1.1
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='{:.1e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    if p<0.05:

        if compr_pos[2] == 't':
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*h,  p_label, ha='center', va='bottom', color=col,fontsize=16) # "{:.1e}".format(p)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*1.1, y2*1.1, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*2, p_label, ha='center', va='bottom', color=col,fontsize=16)

def signal_compr(signals,binding_id,figname,title,colors,cancer_cols,normal_cols):
    
    ctrl_cols = signals.columns.difference(np.append(normal_cols,cancer_cols))
    # normalized binding signal
    cancer_signals = signals.loc[binding_id,cancer_cols].values
    normal_signals = signals.loc[binding_id,normal_cols].values
    ctrl_signals = signals.loc[binding_id,ctrl_cols].values
    
    # plot to compare binding signals
    plt.figure(figsize=(2,2.5))
    box_vals = [ctrl_signals,normal_signals,cancer_signals]
    positions = [1,2,3]
    g = plt.boxplot(box_vals,positions=positions,widths = .55,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='k'),showfliers=False)
    for compr_pos in [[0,2,'t',1.17],[1,2,'t',1.0]]:
        mark_pvalue(compr_pos,positions,box_vals)

    # scatters    
    scatter_X = []
    for position_id in np.arange(len(positions)):
        scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
        plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=20,zorder=0,alpha=1,marker='o')
    
    #plt.axhline(y=0,color='grey',linestyle='--',linewidth=.7)
    plt.ylabel('CTCF binding levels',fontsize=18)
#     plt.axes().tick_params(axis='y',direction='out', length=4, width=.8, colors='black')
#     plt.axes().tick_params(axis='x',direction='out', length=4, width=.8, colors='black')
    #plt.ylim([-5,10])
    #plt.title(title,fontsize=20)
    sns.despine(offset=None, trim=False)
    plt.axes().set_xticklabels(['Others','T cell','T-ALL'],rotation=30,ha='right',fontsize=18) 
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()


   
def main():

    indir='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f5_binding_signal_compr/f1_signal_compr_layout_TALL_cases'
    outdir='figures/cancer_specific_binding_signal_compr_TALL_gained_lost_constitutive'
    os.makedirs(outdir,exist_ok=True)

    signals = pd.read_csv(indir+os.sep+'signals_RPKM_QuantileNormalized.TALL_binding_filtered.csv',index_col=0)
    # ==== cancer datastes and normal datasets for T-ALL ====
    cancercols,normalcols = CTCF_TALL_modules_new.cancer_specific_cancer_normal_gsmID('T-ALL')
    # ==== T-ALL gained and lost bindings
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL');#print(gained_df)
    gained_df = gained_df.sort_values(by=['cancer_vs_other_stats'],ascending=False)
    lost_df = lost_df.sort_values(by=['cancer_vs_other_stats'],ascending=True)
    
    print()
    for binding_id in gained_df.index[:10]:
        title='T-ALL gained'
        colors = [ 'silver','k','red']
        figname = outdir+os.sep+'TALL_gained_signal_compr_{}.pdf'.format(binding_id)
        signal_compr(signals,binding_id,figname,title,colors,cancercols,normalcols)

    for binding_id in lost_df.index[:10]:
        title='T-ALL lost'
        colors = [ 'silver','k','blue']
        figname = outdir+os.sep+'TALL_lost_signal_compr_{}.pdf'.format(binding_id)
        signal_compr(signals,binding_id,figname,title,colors,cancercols,normalcols)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
