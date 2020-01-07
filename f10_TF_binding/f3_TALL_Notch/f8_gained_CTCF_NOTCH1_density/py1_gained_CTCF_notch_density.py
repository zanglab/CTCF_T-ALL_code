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
sns.set(font_scale=1.6)
sns.set_style("whitegrid", {'axes.grid' : False})
import CTCF_TALL_modules_new
import scipy
from scipy import stats
import bisect
sns.set_style("ticks")



def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit');print(s,p)
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),97)*1.05,1.02, 'k'
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



def box_compare_save_file(gained_df_tmp,const_df_tmp,basename):
    
    gained_df_tmp.to_csv(basename+'_gained.csv')
    const_df_tmp.to_csv(basename+'_const.csv')
        
    plt.figure(figsize=(2.3,3))
    a = gained_df_tmp.fillna(0)['tf_num']
    b = const_df_tmp.fillna(0)['tf_num']
    
#     data = [a,b]
    data=[np.log2(a+1),np.log2(b+1)]
    positions=[0,1]
    xticklabels=['gained','Control']
    colors=['red','silver']

#     sns.distplot(b,hist=False,label='Control',color='grey',kde_kws={'cumulative': False})#,kde=False,fit=scipy.stats.gumbel_r,fit_kws={"color": "grey"})#.set(xlim=(0))
#     sns.distplot(a,hist=False,label='T-ALL Lost',color='blue',kde_kws={'cumulative': False})#,kde=False,fit=scipy.stats.gumbel_r,fit_kws={"color": "k"})

    g = plt.boxplot(data,positions=positions,widths = .55,patch_artist=True,\
                 boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='k'),showfliers=False)    
    
#     scatter_X = []
#     for position_id in np.arange(len(positions)):
#         scatter_x = np.random.normal(positions[position_id],0.06,len(data[position_id]))
#         plt.scatter(scatter_x,data[position_id],color=colors[position_id],s=10,zorder=0,alpha=0.99)

    for compr_pos in [[0,1,'t']]:
        mark_pvalue(compr_pos,positions,data)
    #plt.legend(borderaxespad=0.1,labelspacing=.1,fontsize=14)
    
    plt.ylabel('log2(# NOTCH)',fontsize=18)
#     plt.axes().set_ylim([0,8])
    plt.axes().tick_params(axis='x',direction='out', length=4, width=1, colors='black')    
    plt.axes().tick_params(axis='y',direction='out', length=4, width=1, colors='black')  
    sns.despine(offset=0, trim=False)  
    plt.axes().set_xticklabels(xticklabels,rotation=30,fontsize=18,ha='right')
    plt.savefig(basename+'.pdf',bbox_inches='tight',pad_inches=.1,transparent=True,dpi=600)
    plt.close()


def violinplot(gained_df_tmp,const_df_tmp,basename):
    
    plt.figure(figsize=(2.3,3))
    a = gained_df_tmp.fillna(0)['tf_num']
    b = const_df_tmp.fillna(0)['tf_num']
    
#     data = [a,b]
    data=[np.log2(a+1),np.log2(b+1)]
    positions=[0,1]
    xticklabels=['gained','Control']
    colors=['red','silver']

    palette = {0:'red',1:'silver'}
    sns.violinplot(data=data,palette=palette)

    plt.ylabel('log2(# NOTCH)',fontsize=18)
#     plt.axes().set_ylim([0,8])
    plt.axes().tick_params(axis='x',direction='out', length=4, width=1, colors='black')    
    plt.axes().tick_params(axis='y',direction='out', length=4, width=1, colors='black')  
    sns.despine(offset=0, trim=False)  
    plt.axes().set_xticklabels(xticklabels,rotation=30,fontsize=18,ha='right')
    plt.savefig(basename+'_violin.pdf',bbox_inches='tight',pad_inches=.1,transparent=True,dpi=600)
    plt.close()



def count_tf_num(tf_df):

    col_name = tf_df.columns[-1]
    for index in tf_df.dropna().index:
        tf_num = len(tf_df.loc[index,col_name].split(','))
        domain_len = tf_df.loc[index,'domain_100k_1M_right']-tf_df.loc[index,'domain_100k_1M_left']
        tf_df.loc[index,'tf_num'] = int(tf_num)
        tf_df.loc[index,'tf_num_perK'] = 1000*tf_num/domain_len
    return tf_df
    



def main(outdir):

    outdir='f1_density_figs'
    os.makedirs(outdir,exist_ok=True)
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    const_df = CTCF_TALL_modules_new.return_constitutive_df()

#     cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded','PD9','PD30']
    cancertypes=['T-ALL']
    for cancertype in cancertypes: 
        all_feature_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/fu_feature_combination/f2_combined_sep_cancertype/{}_binding_features.csv'.format(cancertype)
        with open(all_feature_file) as all_feature_inf:
            df = pd.read_csv(all_feature_inf,index_col=0)
#         print(df.columns);exit()
        
        # ==== calculate number of intra-domain NOTCH1 binding for gained CTCF and control
        col_names = ['annotation','domain_100k_1M_left', 'domain_100k_1M_right','notch_center']
        gained_tf_density = df.loc[gained_df.index][col_names]
        gained_tf_density = count_tf_num(gained_tf_density)
        const_tf_density = df.loc[const_df.index][col_names]
        const_tf_density = count_tf_num(const_tf_density)
        
        # compare notch binding for all gained and const CTCFs
        gained_df_tmp = gained_tf_density
        const_df_tmp = const_tf_density
        basename = outdir+os.sep+'{}_notch_density_allCTCF'.format(cancertype)
        box_compare_save_file(gained_df_tmp,const_df_tmp,basename)
        violinplot(gained_df_tmp,const_df_tmp,basename)    
        
        # only keep those CTCFs on distal regions
        gained_df_tmp = gained_tf_density.loc[gained_tf_density['annotation']=='Distal']
        const_df_tmp = const_tf_density.loc[const_tf_density['annotation']=='Distal']
        basename = outdir+os.sep+'{}_notch_density_DistalCTCF'.format(cancertype)
        box_compare_save_file(gained_df_tmp,const_df_tmp,basename)
        violinplot(gained_df_tmp,const_df_tmp,basename)
    


    





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
