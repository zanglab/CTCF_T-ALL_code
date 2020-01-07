import os,sys,argparse,glob
import numpy as np
import pandas as pd
import find_overlap_keep_info_NOT_sep_strand_asimport
import re
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=18
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.6)
sns.set_style("whitegrid", {'axes.grid' : False,'grid.color': 'grey'})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})


cancertypes=['BRCA','CRC','LUAD','PRAD','PRAD_TissueAdded']
cancertype_match_names={'BRCA':'BRCA','CRC':'COAD','LUAD':'LUAD','PRAD':'PRAD','PRAD_TissueAdded':'PRAD'}


def compr_signals(df,cancer_columns,ctrl_columns,flag):
    if flag =='raw':
        df = np.log2(df+1)
        
    cancer_df = df[cancer_columns]
    ctrl_df = df[ctrl_columns]
    #s,p = stats.ttest_ind(cancer_df,ctrl_df,axis=1)
    return cancer_df.mean(axis=1)- ctrl_df.mean(axis=1)

def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),100)*0.9 ,1.1, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),0)*0.9
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='{:.1e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    if p<0.05:
#         if p<0.001:
#             p_label = '**'
        if compr_pos[2] == 't':
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*1.25, p_label, ha='center', va='center', color=col,fontsize=14)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*1.1, y2*1.1, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*1.25, p_label, ha='center', va='top', color=col,fontsize=14)
    
def compr_plot_box(gained_df,lost_df,const_df,cancer_columns,ctrl_columns,outdir,cancertype,flag):
    plt.figure(figsize=(2.1,2.6))
    #### double check the function here
    #### changed into z-score
    gained_s = compr_signals(gained_df,cancer_columns,ctrl_columns,flag)
    lost_s = compr_signals(lost_df,cancer_columns,ctrl_columns,flag)
    const_s = compr_signals(const_df,cancer_columns,ctrl_columns,flag)
    box_vals = [lost_s,const_s,gained_s]
    positions = [1,2,3]
    colors = [ 'blue','darkgrey','red']
    
    g = plt.boxplot(box_vals,positions=positions,widths = .55,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='k'),showfliers=False)
    
    scatter_X = []
    for position_id in np.arange(len(positions)):
        scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
        plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=22,zorder=0,alpha=1,marker='o')
    
    for compr_pos in [[0,1,'b'],[1,2,'t']]:
        mark_pvalue(compr_pos,positions,box_vals)

    plt.axhline(y=0,color='grey',linestyle='--',linewidth=.7)
    plt.ylabel('ATAC-seq \ndifferential score'.format(cancertype),fontsize=20)
    #plt.axes().tick_params(axis='x',direction='out', length=0, width=.8, colors='black')
    sns.despine(offset=None, trim=False)
    if cancertype_match_names[cancertype] == 'PRAD':
        plt.ylim([-4,4.9])
    elif cancertype_match_names[cancertype] == 'LUAD':
        plt.ylim([-2,3.5])
    else:
        plt.ylim([-3.5,4.5])

    plt.title('{}'.format(cancertype_match_names[cancertype]),fontsize=20,va='baseline')
    plt.axes().set_xticklabels(['Lost','Control','Gained'],rotation=30,ha='center',fontsize=18) 
    plt.savefig('{}/{}_ATAC_sig_compr.png'.format(outdir,cancertype),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)



def compr_binding_signal(indir,flag,outdir):

    indir = indir+os.sep+flag
    outdir = outdir+os.sep+flag
    os.makedirs(outdir,exist_ok=True)
      
    for cancertype in cancertypes:
        if 1:         
            gained_file = '{}/{}_gained_ATAC_sig.bed'.format(indir,cancertype)
            lost_file = '{}/{}_lost_ATAC_sig.bed'.format(indir,cancertype)
            const_file = '{}/constitutive_ATAC_sig.bed'.format(indir)
                        
            gained_df = pd.read_csv(gained_file,sep='\t')
            gained_df = gained_df.iloc[:,7:]
            lost_df = pd.read_csv(lost_file,sep='\t')
            lost_df = lost_df.iloc[:,7:]
            const_df = pd.read_csv(const_file,sep='\t')
            const_df = const_df.iloc[:,7:]

            cancer_columns = [i for i in gained_df.columns if re.match(cancertype_match_names[cancertype],i)]
            ctrl_columns = [i for i in gained_df.columns if i not in cancer_columns]
            
            compr_plot_box(gained_df,lost_df,const_df,cancer_columns,ctrl_columns,outdir,cancertype,flag)
#             exit()



def main():

    
    indir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f5_TCGA_ATAC/f1_cancer_specific_CTCF_accessibility/f1_ATAC_sig_abstract'
    outdir = 'figures/cancer_specific_TCGA_ATAC_accessibility_compr'
    os.makedirs(outdir,exist_ok=True)
    
    #compr_binding_signal(indir,'raw',outdir)
#     compr_binding_signal(indir,'log2norm',outdir)
    compr_binding_signal(indir,'mynorm',outdir)





           
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--infile1', action = 'store', type = str,dest = 'infile1', help = 'input file to be compared/separated', metavar = '<file>')
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
