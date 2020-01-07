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
sns.set(font_scale=1.4)
sns.set_style("whitegrid", {'axes.grid' : False})
import CTCF_TALL_modules_new
import scipy
from scipy import stats
import bisect
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]



def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit');print(s,p)
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),100)*1.05,1.02, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),0)*1.05
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='{:.1e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    if p<0.05 or np.isnan(p):
        if compr_pos[2] == 't':
            plt.plot([x1+.05, x1+.05, x2-.05, x2-.05], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*1.05, p_label, ha='left', color=col,fontsize=17)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*1.05, y2*1.05, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*1.05, p_label, ha='center', va='top', color=col,fontsize=17)



def box_compare_save_file(a1,b1,c1,df,xticklabels,figname):
    
    a = df.loc[a1]['log2FoldChange'].values
    b = df.loc[b1]['log2FoldChange'].values
    c = df.loc[c1]['log2FoldChange'].values
    
    data=[c,b]
    positions=[1,2]
    colors=['k','red']
#     data=[a,b]
#     positions=[0,1]
#     colors=['silver','red']
    plt.figure(figsize=(2,3))
    g = plt.boxplot(data,positions=positions,widths = .55,patch_artist=True,\
                 boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='k'),showfliers=False)    
    
    scatter_X = []
    for position_id in np.arange(len(positions)):
        scatter_x = np.random.normal(positions[position_id],0.06,len(data[position_id]))
        plt.scatter(scatter_x,data[position_id],color=colors[position_id],s=30,zorder=2,alpha=0.99)

    for compr_pos in [[0,1,'t']]:
        mark_pvalue(compr_pos,positions,data)
    #plt.legend(borderaxespad=0.1,labelspacing=.1,fontsize=14)
    
    plt.ylabel('log$_2$FC',fontsize=18)
#     plt.axes().set_ylim([0,8])
    plt.axes().tick_params(axis='x',direction='out', length=4, width=1, colors='black')    
    plt.axes().tick_params(axis='y',direction='out', length=4, width=1, colors='black')  
#     sns.despine(offset=0, trim=False)  
    plt.axes().set_xticklabels(xticklabels,rotation=30,fontsize=18,ha='right')
    plt.savefig(figname,bbox_inches='tight',pad_inches=.1,transparent=True,dpi=600)
    plt.close()

    

def read_genes(flag,flag2,indir):
    outfile = indir+os.sep+'{}_{}_allgenes.txt'.format(flag2,flag)
    with open(outfile,'r') as outf:
        genes = [ii.strip() for ii in outf.readlines()]
    
    outfile = indir+os.sep+'{}_{}_shCTCF_dngenes.txt'.format(flag2,flag)
    with open(outfile,'r') as outf:
        genes2 = [ii.strip() for ii in outf.readlines()]
        
    return genes, genes2


def main(outdir):

    indir = 'f3_dNOTCH_gainedCTCF_expr_in_shCTCF'
    outdir = 'f4_dNOTCH_gainedCTCF_expr_in_shCTCF_figs'
    os.makedirs(outdir,exist_ok=True)

    deseq2_file = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/18_MYC-ChIP_shCTCF-RNA/shCTCF_RNA/f0_processing/salmon_Deseq2/salmon_Deseq2_pca/CUTLL1/f4_deseq_out_shrink/treated_CUTLL1_shCTCF_vs_ctrl_CUTLL1_PIG.csv"
    with open(deseq2_file) as inf:
        df = pd.read_csv(inf,index_col=0)
    
    flags = ['all_genes','dNotch_hic_increased_gained_ctcf_intra_domain','dNotch_hic_NOT_increased_gained_ctcf_intra_domain']
    xticklabels = ['all','dNOTCH w/ gained CTCF\n(increased Hi-C contact)']
    xticklabels = ['De-interaction','In-interaction']
    
    flag2 = 'cutll1_upgenes'
    a1,a2 = read_genes(flags[0],flag2,indir)#;print(len(a1),len(a2))
    b1,b2 = read_genes(flags[1],flag2,indir)
    c1,c2 = read_genes(flags[2],flag2,indir)
    
    figname = outdir+os.sep+'dNOTCH_gained_CTCF_hic_increased_shCTCF_expr.pdf'
    box_compare_save_file(a1,b1,c1,df,xticklabels,figname)
    





if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
#     parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.outdir)
