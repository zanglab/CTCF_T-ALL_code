import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import association_with_genes
import association_with_regions
import re,bisect
import CTCF_TALL_modules_new
import scipy
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=15
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.5)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")

    
def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit');print(s,p)
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),99)*1.05,1.02, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),2)*0.35
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='{:.1e}'.format(p)
    if p<0.05 or np.isnan(p):
        if compr_pos[2] == 't':
            plt.plot([x1+.05, x1+.05, x2-.05, x2-.05], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*1.05, p_label, ha='center', color=col,fontsize=17)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*.91, y2*.91, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*.85, p_label, ha='center', va='top', color=col,fontsize=17)

def compr_gene_len(master_df,figname,dis_marker,distal_thre=3):  
    '''
    compare length of gene body
    ''' 
    fig = plt.figure(figsize=(3.2,3))
    # gene length
    proximal_notch_gene_len = master_df.loc[master_df['{}_distal_type'.format(dis_marker)]=='Notch_Proximal','gene_len']
    distal_notch_distal_CTCF_gene_len = master_df.loc[master_df['{}_distal_type'.format(dis_marker)]=='Notch_Distal_CTCF_Distal','gene_len']
    distal_notch_proximal_CTCF_gene_len = master_df.loc[master_df['{}_distal_type'.format(dis_marker)]=='Notch_Distal_CTCF_Proximal','gene_len']
    
    data=[proximal_notch_gene_len,distal_notch_distal_CTCF_gene_len,distal_notch_proximal_CTCF_gene_len]
    data=[np.log10(proximal_notch_gene_len),np.log10(distal_notch_distal_CTCF_gene_len),np.log10(distal_notch_proximal_CTCF_gene_len)]

    positions=[0,1,2]
    xticklabels=['Proximal Notch','Distal Notch\nDistal CTCF','Distal NOTCH\nProximal CTCF']
    colors=['b','r','g']
    g = plt.boxplot(data,positions=positions,widths = .55,patch_artist=True,\
                 boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='k'),showfliers=False)    
    
    scatter_X = []
    for position_id in np.arange(len(positions)):
        scatter_x = np.random.normal(positions[position_id],0.06,len(data[position_id]))
        plt.scatter(scatter_x,data[position_id],color=colors[position_id],s=30,zorder=0,alpha=0.99)

    for compr_pos in [[0,1,'t']]:
        mark_pvalue(compr_pos,positions,data)
    
    plt.ylabel('Length of GB (log10 bps)',fontsize=18)
#     plt.axes().set_ylim([0,8])
    plt.axes().tick_params(axis='x',direction='out', length=4, width=1, colors='black')    
    plt.axes().tick_params(axis='y',direction='out', length=4, width=1, colors='black')  
    sns.despine(offset=0, trim=False)  
    plt.axes().set_xticklabels(xticklabels,rotation=45,fontsize=16,ha='right')
    plt.savefig(figname,bbox_inches='tight',pad_inches=.1,transparent=True,dpi=600)
    plt.close()



def ctcf_notch_gene_dis(master_df,figname,dis_marker,distal_thre=3,labelMark=False):
    ''' for each CTCF-notch-gene pair, plot distance between 
    ctcf-gene against notch-gene
    and assign distal/proximal type for each index/pair
    ''' 
    fig = plt.figure(figsize=(4,4))
    for ii in master_df.index:
#         notch_dis = master_df.loc[ii,'dis_notch_gene']
#         ctcf_dis = master_df.loc[ii,'dis_ctcf_gene']
        notch_dis = np.log10(master_df.loc[ii,'dis_notch_{}'.format(dis_marker)]+1)
        ctcf_dis = np.log10(master_df.loc[ii,'dis_ctcf_{}'.format(dis_marker)]+1)
        gene_name = master_df.loc[ii,'gene_id']
        if notch_dis<distal_thre:
            master_df.loc[ii,'{}_distal_type'.format(dis_marker)]='Notch_Proximal'
            aa = plt.scatter(notch_dis,ctcf_dis,color='b',s=50)
        elif ctcf_dis>distal_thre:
            master_df.loc[ii,'{}_distal_type'.format(dis_marker)]='Notch_Distal_CTCF_Distal'
            bb = plt.scatter(notch_dis,ctcf_dis,color='r',s=50)
        else:
            master_df.loc[ii,'{}_distal_type'.format(dis_marker)]='Notch_Distal_CTCF_Proximal'
            cc = plt.scatter(notch_dis,ctcf_dis,color='g',s=50)
        if labelMark:
            plt.axes().text(notch_dis,ctcf_dis,gene_name,fontsize=6)    
    plt.legend([aa,bb,cc],['Proximal NOTCH','Distal NOTCH\nDistal CTCF','Distal NOTCH\nProximal CTCF'],bbox_to_anchor=[1,1],fontsize=16,frameon=False,markerscale=1.5,\
               borderaxespad=1,labelspacing=1,handletextpad=0.2,handlelength=1,ncol=1)
    sns.despine(offset=0, trim=False)
    plt.xlabel('NOTCH gene distance\n(log10 bps)',fontsize=18)
    plt.ylabel('CTCF gene distance\n(log10 bps)',fontsize=18)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()
    return master_df

    


def main():

    outdir = 'f7_ctcf_notch_gene_figs'
    os.makedirs(outdir,exist_ok=True) 

    master_df = pd.read_csv('f6_ctcf_notch_gene_master/T_ALL_gained_CTCF_NOTCH_gene_master.csv',index_col=0)
    distal_thre=np.log10(2000)
    
    dis_markers = ['gene','promoter']
    for dis_marker in dis_markers:
        # plot distance of notch/ctcf to gene
        figname = outdir+os.sep+'ctcf_notch_{}_dis.pdf'.format(dis_marker)
        ctcf_notch_gene_dis(master_df,figname,dis_marker,distal_thre=distal_thre)
        figname = outdir+os.sep+'ctcf_notch_{}_dis_withMark.pdf'.format(dis_marker)
        master_df= ctcf_notch_gene_dis(master_df,figname,dis_marker,distal_thre=distal_thre,labelMark=True)
        
        figname = outdir+os.sep+'gb_len_by_distal_proximal_notch_{}_dis.pdf'.format(dis_marker)
        compr_gene_len(master_df,figname,dis_marker,distal_thre=distal_thre)
    
    master_df.to_csv(outdir+os.sep+'T_ALL_gained_CTCF_NOTCH_gene_master_new.csv')
        
          
            


 
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

