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
matplotlib.rcParams["font.sans-serif"] = ["Arial"]


def dis_compare(gained_df,lost_df,const_df,df,figname,flag):
    plt.figure(figsize=(4.6,3.2))
    sns.distplot(df.loc[const_df.index]['nearest_pos'],hist=False,label='Control',color='grey',kde_kws={'cumulative': False})#,kde=False,fit=scipy.stats.gumbel_r,fit_kws={"color": "grey"})#.set(xlim=(0))
    sns.distplot(df.loc[lost_df.index]['nearest_pos'],hist=False,label='T-ALL Lost',color='blue',kde_kws={'cumulative': False})#,kde=False,fit=scipy.stats.gumbel_r,fit_kws={"color": "k"})
    sns.distplot(df.loc[gained_df.index]['nearest_pos'],hist=False,label='T-ALL Gained',color='red',kde_kws={'cumulative': False})#,kde=False,fit=scipy.stats.gumbel_r,fit_kws={"color": "firebrick"})
    plt.legend(bbox_to_anchor=[.8,1.1],loc='upper left',borderaxespad=0.3,labelspacing=.3,fontsize=17,handlelength=1.2,handletextpad=0.3,frameon=False)
    
    plt.axes().tick_params(axis='x',direction='out', length=4, width=1, colors='black')    
    plt.axes().tick_params(axis='y',direction='out', length=4, width=1, colors='black')  
    plt.xlabel('Distance to nearest {} (bp)'.format(flag),fontsize=20)  
    plt.ylabel('PDF',fontsize=20)
    plt.text(8.5,-0.06,'log$_{{10}}$',fontsize=19)
    sns.despine(offset=0, trim=False)
    plt.axes().set_xlim([0,8])
    plt.savefig(figname,bbox_inches='tight',pad_inches=.1,transparent=True,dpi=600)
    plt.close()


def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit');print(s,p)
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),99.9)*1.05,1.02, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),1)*0.35
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='*'
    if p<0.05 or np.isnan(p):
        if p<0.001 or np.isnan(p):
            p_label = '**'
        if compr_pos[2] == 't':
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*1.01, p_label, ha='center', color=col,fontsize=25)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*.91, y2*.91, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*.85, p_label, ha='center', va='top', color=col,fontsize=25)


def box_compare(gained_df,lost_df,const_df,df,figname,flag):
    
    plt.figure(figsize=(2.3,3))
    a = df.loc[const_df.index]['nearest_pos']
    b = df.loc[lost_df.index]['nearest_pos']
    c = df.loc[gained_df.index]['nearest_pos']
    data=[b,a,c]
    positions=[0,1,2]
    xticklabels=['T-ALL lost','Control','T-ALL gained']
    colors=['blue','silver','red']
    
#     g = sns.swarmplot(data=data,size=1.7,palette=colors,linewidth=.1)
#     g = sns.boxplot(data=data,width=0.55,linewidth=1.6)
#     for i,box in enumerate(g.artists):
#        box.set_edgecolor('k')
#        box.set_facecolor('w')

    g = plt.boxplot(data,positions=positions,widths = .55,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='k'))    
    
    colors = [ 'blue','silver','red']
    scatter_X = []
    for position_id in np.arange(len(positions)):
        scatter_x = np.random.normal(positions[position_id],0.07,len(data[position_id]))
        plt.scatter(scatter_x,data[position_id],color=colors[position_id],s=20,zorder=0,alpha=0.99)

    
#     g = plt.boxplot(data,positions=positions,widths=0.55,patch_artist=True)
#     for patch, color in zip(g['boxes'], colors):
#         patch.set_facecolor('w')
#         patch.set_edgecolor('k')
    
    for compr_pos in [[0,1,'t'],[1,2,'t'],[0,2,'t']]:
        mark_pvalue(compr_pos,positions,data)
    #plt.legend(borderaxespad=0.1,labelspacing=.1,fontsize=14)
    plt.text(-.9,8.5,'$log10$',fontsize=16)
    plt.ylabel('Distance to \nnearest {}'.format(flag),fontsize=20)
    plt.axes().set_ylim([0,8])
    plt.axes().tick_params(axis='x',direction='out', length=4, width=1, colors='black')    
    plt.axes().tick_params(axis='y',direction='out', length=4, width=1, colors='black')  
    sns.despine(offset=0, trim=True)  
    plt.axes().set_xticklabels(xticklabels,rotation=30,fontsize=18,ha='right')
    plt.savefig(figname,bbox_inches='tight',pad_inches=.1,transparent=True,dpi=600)
    plt.close()


def main():

    outdir='f2_dis_figs'
    os.makedirs(outdir,exist_ok=True)
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    const_df = CTCF_TALL_modules_new.return_constitutive_df()

    flags = {"union_CTCF_nearest_CUTLL1_MYC":"MYC","union_CTCF_nearest_NOTCH1_dynamic":"NOTCH1-D","union_CTCF_nearest_NOTCH1_w4h":"NOTCH1"}
    infiles = glob.glob('f1_NOTCH1_MYC_nearest_dis_csv/*.csv')
    for infile in infiles:
        basename = os.path.basename(infile).split('.csv')[0];print(basename)
        with open(infile) as inf:
            df = pd.read_csv(inf,index_col=3,sep='\t')
        df['nearest_pos']= pd.concat([-df['left_pos']+df['start'],df['right_pos']-df['start']],axis=1).min(axis=1)
        df['nearest_pos'] = np.log10(df['nearest_pos'])
        
        figname='{}/{}_dis.pdf'.format(outdir,basename)
        dis_compare(gained_df,lost_df,const_df,df,figname,flags[basename])

        figname='{}/{}_box.png'.format(outdir,basename)
       # box_compare(gained_df,lost_df,const_df,df,figname,flags[basename])#;exit()





if __name__ == '__main__':

  
    main()
