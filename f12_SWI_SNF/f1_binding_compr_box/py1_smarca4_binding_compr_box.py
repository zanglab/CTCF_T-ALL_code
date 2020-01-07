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
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),99.9)*1.1 ,1.1, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),0)*1
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='{:.1e}'.format(p)
#     if p_label[-2]=='0':
#         p_label = p_label[:-2]+p_label[-1]
    if p<0.05:

        if compr_pos[2] == 't':
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*h,  p_label, ha='center', va='bottom', color=col,fontsize=16) # "{:.1e}".format(p)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2-0.2, y2-0.5, y2-0.5, y2-0.2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2-1.8, p_label, ha='center', va='bottom', color=col,fontsize=16)

def signal_compr(rpkm_df,gained_df,lost_df,const_df,basename,title,colors,cancertype):
    
    # normalized binding signal
    gained_signals = rpkm_df.loc[gained_df.index,1];gained_signals.to_csv(basename+'_gained.csv')
    lost_signals = rpkm_df.loc[lost_df.index,1];lost_signals.to_csv(basename+'_lost.csv')
    ctrl_signals = rpkm_df.loc[const_df.index,1];ctrl_signals.to_csv(basename+'_ctrl.csv')
    
    # plot to compare binding signals
    plt.figure(figsize=(2,2.5))
    box_vals = [lost_signals.values,ctrl_signals.values,gained_signals.values]
    positions = [1,2,3]
    g = plt.boxplot(box_vals,positions=positions,widths = .55,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='k'),showfliers=False)
    for compr_pos in [[0,1,'b',1.1],[1,2,'t',1.0]]:
        mark_pvalue(compr_pos,positions,box_vals)

    # scatters    
    scatter_X = []
    for position_id in np.arange(len(positions)):
        scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
        plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=20,zorder=0,alpha=1,marker='o')
    
    #plt.axhline(y=0,color='grey',linestyle='--',linewidth=.7)
    plt.ylabel('SMARCA4\nbinding levels',fontsize=18)
#     plt.axes().tick_params(axis='y',direction='out', length=4, width=.8, colors='black')
#     plt.axes().tick_params(axis='x',direction='out', length=4, width=.8, colors='black')
    plt.ylim([-2,8])
    #plt.title(title,fontsize=20)
    sns.despine(offset=None, trim=False)
    plt.axes().set_xticklabels(['{} Lost'.format(cancertype),'Constitutive','{} Gained'.format(cancertype)],rotation=30,ha='right',fontsize=18) 
    plt.savefig(basename+'_box.pdf',bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()


   
def main():

    outdir='f1_figs'
    os.makedirs(outdir,exist_ok=True)
    rpkm_file_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f12_SWI_SNF/f0_SMARCA4_union_binding/binding_RPKM/"
    colors = [ 'blue','silver','red']
    
    # ==== AML case
    cancertype = 'AML'
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype);#print(gained_df)
    const_df = CTCF_TALL_modules_new.return_constitutive_df()
    # ==== RPKM binding_region
    rpkm_file = rpkm_file_dir+os.sep+'EOL1_SMARCA4_binding_region_RPKM.csv'
    rpkm_df = pd.read_csv(rpkm_file,index_col=0,sep='\t',header=None)
    rpkm_df = np.sqrt(rpkm_df)
    title='AML gained'
    basename = outdir+os.sep+'AML_EOL1_SMARCA4_binding'
    signal_compr(rpkm_df,gained_df,lost_df,const_df,basename,title,colors,cancertype)


#     cancertype = 'CRC'
#     gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype);#print(gained_df)
#     const_df = CTCF_TALL_modules_new.return_constitutive_df()
#     # ==== RPKM binding_region
#     rpkm_file = rpkm_file_dir+os.sep+'HCT116_SMARCA4_binding_region_RPKM.csv'
#     rpkm_df = pd.read_csv(rpkm_file,index_col=0,sep='\t',header=None)
#     rpkm_df = np.sqrt(rpkm_df)
#     title='AML gained'
#     basename = outdir+os.sep+'CRC_HCT116_SMARCA4_binding'
#     signal_compr(rpkm_df,gained_df,lost_df,const_df,basename,title,colors,cancertype)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
