import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
import association_with_genes
import association_with_regions
import re,bisect
import CTCF_TALL_modules_new
import scipy
from scipy import stats
sys.path.insert(0,os.path.abspath('modules'))
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"


def plot_dynamic_scatter(df,df_tmp,columns,outdir,basename,xlabel,ylabel,title):

    fig = plt.figure(figsize=(2.8,2.8))
    figname = outdir+os.sep+'{}.pdf'.format(basename)
    g1 = plt.scatter(df[columns[0]],df[columns[1]],c='silver',s=30,alpha=1)
    g2 = plt.scatter(df_tmp[columns[0]],df_tmp[columns[1]],c='red',s=30,alpha=1)
    plt.axhline(y=0,color='gray',linestyle='--',linewidth=.7)
    plt.axis('equal')
#     sns.despine(offset=0, trim=False)
    
    plt.xlabel(xlabel,fontsize=15)
    plt.ylabel(ylabel,fontsize=15)
    plt.xlim([-1.5,4.2])   
    plt.ylim([-2,2])
    legend = plt.legend([g1,g2],['T-ALL CTCF sites','T-ALL gained'],loc='upper left',bbox_to_anchor=[-.1,1.25],markerscale=1.5,fontsize=14,borderaxespad=0.1,labelspacing=.2,handletextpad=0.1,frameon=False)
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_edgecolor("white")
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,dpi=600,transparent=True)
    plt.close()


def mark_pvalue(compr_pos,positions,box_vals,flag):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    print('\n',flag,s,p)
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),99)*1.1 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),0)*0.95
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]

    p_label='{:.1e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    if p<0.05:
        if p<0.05:
            star_mark="*"
            if p<0.001:
                star_mark="**"

        if compr_pos[2] == 't':
            plt.plot([x1, x1, x2, x2], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*1.4, star_mark, ha='center', va='bottom', color=col,fontsize=18)
        else:
            plt.plot([x1, x1, x2, x2], [y2, y2*1.1, y2*1.1, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*1.3, star_mark, ha='center', va='bottom', color=col,fontsize=18)

def box_plot_for_scatter(df,df_tmp,columns,outdir,basename,xlabel,ylabel):

    fig = plt.figure(figsize=(1.5,2.8))
    figname = outdir+os.sep+'{}_box.pdf'.format(basename)
    positions = [1,2]    
    plot_lists = [df[columns[1]],df_tmp[columns[1]]]
    g = plt.boxplot(plot_lists,positions=positions,widths = .5,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)    
    for compr_pos in [[0,1,'t']]:
        mark_pvalue(compr_pos,positions,plot_lists,basename)
    colors = [ 'silver','red']
    scatter_X = []
    for position_id in np.arange(len(positions)):
        scatter_x = np.random.normal(positions[position_id],0.07,len(plot_lists[position_id]))
        plt.scatter(scatter_x,plot_lists[position_id],color=colors[position_id],s=20,zorder=0,alpha=0.99)

#     sns.despine(offset=0, trim=False)
    plt.axes().set_xticks(positions)
    #plt.legend([g["boxes"][0],g["boxes"][1]],['All','Gained'],loc='upper right',fontsize=12)
    plt.ylabel(ylabel,fontsize=15)
    plt.ylim([-2.5,2.5])
    
    plt.axes().set_xticklabels(['T-ALL CTCF sites','T-ALL gained'],rotation=30,ha='right',fontsize=15) 
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,dpi=600,transparent=True)
    plt.close()


def notch_overlap_info_by_group(df,df_tmp,outdir,basename):

    columns = ['gsi_over_dmso_log_avg', 'gsi_over_dmso_adj_logFC']
    xlabel = 'Average log$_2$(RPKM)'
    ylabel = 'CTCF ChIP-seq log$_2$(fold change)\nGSI-3d over DMSO'
    title="GSI vs. DMSO"
    plot_dynamic_scatter(df,df_tmp,columns,outdir,basename+'_gsi_over_dmso',xlabel,ylabel,title)
    box_plot_for_scatter(df,df_tmp,columns,outdir,basename+'_gsi_over_dmso',xlabel,ylabel)
    
    columns = ['w4hr_over_gsi_log_avg', 'w4hr_over_gsi_adj_logFC']
    xlabel = 'Average log$_2$(RPKM)'
    ylabel = 'CTCF ChIP-seq log$_2$(fold change)\nGSI-washout over GSI'
    title="GSI washout"
    plot_dynamic_scatter(df,df_tmp,columns,outdir,basename+'_w4hr_over_gsi',xlabel,ylabel,title)
    box_plot_for_scatter(df,df_tmp,columns,outdir,basename+'_w4hr_over_gsi',xlabel,ylabel)
    


def main():

    outdir = 'f3_tall_gained_GSI_bidings_changes'
    os.makedirs(outdir,exist_ok=True) 
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    
    infile = 'f2_dynamic_MA_csv_log2/Jurkat_GSI_union_bindings_MA_adj_logFC.csv'
    with open(infile) as inf:
        df = pd.read_csv(inf,index_col =0)    

    df_tmp = df.loc[df.index.intersection(gained_df.index)]
    
    basename = 'gsi_peak_overlapped'
    notch_overlap_info_by_group(df,df_tmp,outdir,basename)




 
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

