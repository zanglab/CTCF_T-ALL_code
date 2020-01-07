import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=15
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
import association_with_genes
import association_with_regions
import re,bisect
# import CTCF_TALL_modules
import scipy
from scipy import stats
# import boxplot_with_pvalue

sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]

def mark_pvalue(compr_pos,positions,box_vals,flag):

    s,p = stats.ttest_rel(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    print('\n',flag,s,p)
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),98)*1.1 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),0)*0.95
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    #if compr_pos[0]%3==0:
    #    y = y*1.5
    p_label='{:.1e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    p_label='*'
    if p<0.001:
        p_label = '**'
    
    if p<0.05:

        if compr_pos[2] == 't':
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*1.05, p_label, ha='center', color=col,fontsize=15)
        else:
            plt.plot([x1, x1, x2, x2], [y2, y2*1.1, y2*1.1, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*1.3,p_label, ha='center', va='bottom', color=col,fontsize=14)



def box_plot_for_zscores(plot_lists,figname,xticklabels=None,legends=None,flag=None):

    fig = plt.figure(figsize=(2,2.7))
    positions = [1,2]    
    g = plt.boxplot(plot_lists,positions=positions,widths = .52,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)    

    for compr_pos in [[0,1,'t']]:
        mark_pvalue(compr_pos,positions,plot_lists,flag)
    
    colors = [ 'darkgrey','red']
    scatter_X = []
    for position_id in np.arange(len(positions)):
        scatter_x = np.random.normal(positions[position_id],0.07,len(plot_lists[position_id]))
        plt.scatter(scatter_x,plot_lists[position_id],color=colors[position_id],s=30,zorder=0,alpha=0.99)

    plt.axes().set_xticks(np.array(positions))
    #plt.axhline(y=0,color='grey',linestyle='--',linewidth=.7)
    #plt.axvline(x=0,color='grey',linestyle='--',linewidth=.5)
    #plt.legend([g["boxes"][0],g["boxes"][1]],legends,loc='upper right',fontsize=12)
    #plt.ylabel(columns[1])
    #if re.search('lost',flag):
    plt.ylabel('Hi-C (z-score) of\nCTCF and NOTCH1',fontsize=18)
    plt.title((''.join(flag.split('_')[-1]).capitalize()),fontsize=18)
    #plt.axes().tick_params(axis='y',direction='out', length=3, width=.7, colors='black')
    plt.ylim([-2,3])
#     plt.axes().set_yticks([])
#     sns.despine(offset=0, trim=False)
    plt.axes().set_xticklabels(xticklabels,rotation=0,ha='center',fontsize=18) 
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,dpi=600,transparent=True)
    plt.close()


def return_notch_bg_scores(score_file):

    with open(score_file) as score_inf:
        score_df = pd.read_csv(score_inf,sep='\t')
    #print(score_df);exit()       
    compr_zscore = []
    for id_index in score_df.index:
        compr_score = score_df.loc[id_index,'ctcf_notch_interaction']
        bg_score = [float(i) for i in score_df.loc[id_index,'intra_domain_ctcf_all_interactions'].split(',')]
        z_scores = stats.zscore(np.append(compr_score,bg_score))
        compr_zscore.append(z_scores[0])

    return score_df,compr_zscore




def main():

    outdir = 'f4_interaction_compr_zscores'
    os.makedirs(outdir,exist_ok=True) 

    for basename in ['T-ALL_gained','T-ALL_lost']:
        for hic_normalization in ['raw']:
            a6010_score_file = 'f2_CTCF_Notch1_hic_interaction/{}_ctcf_notch_interaction_vs_bg_in_a6010_{}.csv'.format(basename,hic_normalization)
            cutll1_score_file = 'f2_CTCF_Notch1_hic_interaction/{}_ctcf_notch_interaction_vs_bg_in_cutll1_{}.csv'.format(basename,hic_normalization)          
            print(basename,hic_normalization,os.path.isfile(a6010_score_file),os.path.isfile(cutll1_score_file))
            a_df,a_zscores = return_notch_bg_scores(a6010_score_file)#;print(a_df);exit()
            j_df,j_zscores = return_notch_bg_scores(cutll1_score_file);print(a_df,j_df)
            
            # print the z-score change for each notch
            outfile=outdir+os.sep+'{}_{}_zscore.txt'.format(hic_normalization,basename)
            outf = open(outfile,'w')
            outf.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('id','chr','ctcf_pos','notch_pos','cd4_zscore','cutll1_zscore','delta_zscore'))
            
            for index in a_df.index:
                id = a_df.loc[index,'id']
                chr = a_df.loc[index,'chr']
                ctcf_pos = a_df.loc[index,'ctcf_pos']
                notch_pos = a_df.loc[index,'notch_pos']
                a_zscore = a_zscores[index]
                j_zscore = j_zscores[index]
                outf.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(id,chr,ctcf_pos,notch_pos,a_zscore,j_zscore,j_zscore-a_zscore))
            outf.close()    
    
            z_score_df = pd.read_csv(outfile,sep='\t')
            z_score_df = z_score_df.sort_values(by='delta_zscore',ascending=False).drop_duplicates('id',keep='first')
            z_score_df.to_csv(outdir+os.sep+'{}_{}_zscore_delta_change_filtered.txt'.format(hic_normalization,basename),index=False,sep='\t')
            ## plot the z-score figs
            a_zscores,j_zscores = z_score_df['cd4_zscore'],z_score_df['cutll1_zscore']
            s,p = stats.ttest_rel(a_zscores,j_zscores)
            #print(a_zscores,j_zscores);exit()
            print(basename,hic_normalization,s,p)
            figname = outdir+os.sep+'{}_{}_zscore.pdf'.format(hic_normalization,basename)
            plot_lists = [a_zscores,j_zscores]
            xticklabels = ['T-cell','T-ALL']
            legends = ['ctcf w/ intra-domain','ctcf w/ notch']
            flag = '{}'.format(basename,hic_normalization)
            box_plot_for_zscores(plot_lists,figname,xticklabels = xticklabels,legends = legends,flag = flag)




 
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

