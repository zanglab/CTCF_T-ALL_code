import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
#import myplot
import CTCF_TALL_modules_new
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"


chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']

cancer_type_hic_filename = {'T-ALL':['Jurkat','A6010'],\
                                'T-ALL-P1':['PD9','A6010'],\
                                'T-ALL-P2':['PD31','A6010'],\
                                'CRC':['HCT116','trans_colon1'],\
                                'PRAD':['LNCaP','PrEC'],\
                                'PRAD_TissueAdded':['LNCaP','PrEC']}

cancertype_title = {'T-ALL':'T-ALL','T-ALL-P1':'T-ALL patient1','T-ALL-P2':'T-ALL patient2',\
                    'CRC':'CRC','PRAD':'PRAD','PRAD_TissueAdded':'PRAD'}

pvalue_type,stats_type = 'pvalue','stats'
#pvalue_type,stats_type = 'log_pvalue','log_stats'

################
# for scatter plot
################
def plot_each_point(cancertype,plot_df,color):
    #plot_df = plot_df.dropna();#print(plot_df);exit()
    print(cancertype,plot_df.shape)
    for index in plot_df.index:
        binding = plot_df.loc[index]
        cancer_mean = '{}_mean'.format(cancer_type_hic_filename[cancertype][0])
        normal_mean = '{}_mean'.format(cancer_type_hic_filename[cancertype][1])
        g = plt.scatter(np.log2(binding[cancer_mean]/binding[normal_mean]),-1*np.log10(binding[pvalue_type]),color = color,s=15,marker='o')
    return g
    
def scatter_plot(gained_df,lost_df,const_df,cancertype,normalization,resolution,viewregion,outdir):
    
    figname = outdir+os.sep+'{}_{}_res{}_view{}.pdf'.format(cancertype,normalization,resolution,viewregion)
    title=cancertype_title[cancertype]

    fig = plt.figure(figsize=(3,3))
    a = plot_each_point(cancertype,const_df,'silver')    
    b = plot_each_point(cancertype,gained_df,'red') 
    c = plot_each_point(cancertype,lost_df,'blue') 
    legend = plt.legend([a,c,b,],['Constitutive','{} lost'.format(title.split(' ')[0]),'{} gained'.format(title.split(' ')[0])],\
             markerscale=1.5,fontsize=14,loc='upper left',bbox_to_anchor=[-0.05,1],frameon=False,borderaxespad=0 ,labelspacing=.0,handletextpad=0.1)
    legend.get_frame().set_facecolor('w')
    legend.get_frame().set_edgecolor('w')
    plt.axhline(y=-1*np.log10(0.05),color='grey',linestyle='--',linewidth=.7)
    plt.ylabel('-log$_{{10}}$ ($P$ value)',fontsize=18)
#     plt.xlabel('log$_2$ FC',fontsize=20)
    plt.xlabel('$\Delta$(chromatin interaction)\n{} over normal'.format(title),fontsize=18)
#     plt.title(title,fontsize=18)
    plt.xlim([-3,3])
    plt.axes().tick_params(zorder=20)
    plt.rcParams["axes.axisbelow"] = False
#     sns.despine(offset=None, trim=False)
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .1,dpi=600,transparent=True)
    plt.close()


################
# for box plot
################
def return_sig_nums(cancertype,df,pvalue_thre=0.05):

    total = len(df.index)
    pos_num = len(df.loc[(df[stats_type]>0)&(df[pvalue_type]<pvalue_thre)].index)
    minus_num = len(df.loc[(df[stats_type]<0)&(df[pvalue_type]<pvalue_thre)].index)
    
    df = df.loc[df[pvalue_type]<0.05]
    cancer_mean = '{}_mean'.format(cancer_type_hic_filename[cancertype][0])
    normal_mean = '{}_mean'.format(cancer_type_hic_filename[cancertype][1])
    values = np.log2(df[cancer_mean]/df[normal_mean])
    values = values.replace([np.inf,-np.inf],np.nan).dropna(0)
    #print(values)
    #values = [value for value in values if (value!=np.inf and value!=-np.inf)]
    
    return total,pos_num,minus_num,values

def mark_pvalue(compr_pos,positions,box_vals,figname):

    s,p = stats.ttest_ind(box_vals[compr_pos[1]],box_vals[compr_pos[0]],nan_policy='omit')
    print('\n',os.path.basename(figname),s,p)
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),99)*1.05 ,1.15, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),1)*1.05
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    #if compr_pos[0]%3==0:
    p_label='{:.1e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    p_label='*'
    if p<0.001:
        p_label = '**'
    
    if p<0.05:
        if compr_pos[2] == 't':
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*h, p_label, va='bottom', color=col,fontsize=17)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*1.1, y2*1.1, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*1.2-0.5 , p_label, ha='center',  va='bottom', color=col,fontsize=17)

def box_plot(cancertype,gained_values,lost_values,const_values,normalization,resolution,viewregion,figname):
    
    fig = plt.figure(figsize=(3,3))
    box_vals = [lost_values,const_values,gained_values]
    positions = [1,2,3]
    g = plt.boxplot(box_vals,positions = positions,widths = .6,patch_artist=True,\
        #boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
        medianprops=dict(color='k'),showfliers=False)
    #s,p = stats.ttest_ind(gained_values,const_values);print(s,p)
    colors = [ 'blue', 'darkgrey','red']
    for patch, color in zip(g['boxes'], colors):
        patch.set_facecolor(color)
    
    for compr_pos in [[0,1,'b'],[1,2,'t']]:
        mark_pvalue(compr_pos,positions,box_vals,figname)
    
#     scatter_X = []
#     for position_id in np.arange(len(positions)):
#         scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
#        # plt.scatter(scatter_x,box_vals[position_id],color='k',s=14,zorder=10,alpha=0.9,marker='o',facecolors='none')
#         plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=14,zorder=10,alpha=0.9,marker='o')

    plt.axhline(y=0,color='grey',linestyle='--',linewidth=.7)
#     plt.ylabel('log$_2$ FC',fontsize=20)  
#     plt.xlim([0.0,4])
#     if cancertype=='T-ALL':
    plt.ylim([-2.6,2.6])
    title = cancertype_title[cancertype]
    plt.ylabel('$\Delta$(chromatin interaction)\n{} over normal'.format(title),fontsize=18)
#     plt.title(title,fontsize=18)
    #plt.axes().tick_params(axis='y',direction='out', length=3, width=.8, colors='black')
#     sns.despine(offset=None, trim=False)
    #plt.axes().set_xticklabels(['{}$_{{lost}}$'.format(title.split(' ')[0]),'Constitutive','{}$_{{gained}}$'.format(title.split(' ')[0])],rotation=30,ha='center',fontsize=18) 
    plt.axes().set_xticklabels(['{} lost'.format(title.split(' ')[0]),'Constitutive','{} gained'.format(title.split(' ')[0])],rotation=30,ha='center',fontsize=18) 
    plt.savefig(figname,bbox_inches = 'tight',pad_inches = .05,dpi=600,transparent=True)
    plt.close()



def read_ttest_df(indir,cancertype,normalization,resolution,viewregion):

    gained_df = pd.read_csv(indir+os.sep+'{}_{}_{}_res-{}_viewregion-{}_paired_ttest.csv'.format(cancertype,'gained',normalization,resolution,viewregion),index_col=0)
    gained_df = gained_df.fillna(0)
    lost_df = pd.read_csv(indir+os.sep+'{}_{}_{}_res-{}_viewregion-{}_paired_ttest.csv'.format(cancertype,'lost',normalization,resolution,viewregion),index_col=0)
    lost_df = lost_df.fillna(0)
    const_df = pd.read_csv(indir+os.sep+'{}_{}_{}_res-{}_viewregion-{}_paired_ttest.csv'.format(cancertype,'ctrl',normalization,resolution,viewregion),index_col=0)
    const_df = const_df.fillna(0)

    if cancertype in ['T-ALL-P1','T-ALL-P2']:
        gained_df_new,lost_df_new = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    else:
        gained_df_new,lost_df_new = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype)
        
    gained_df = gained_df.loc[gained_df_new.index]
    lost_df = lost_df.loc[lost_df_new.index]
    print(cancertype,'gained',gained_df.shape)
    print(cancertype,'lost',lost_df.shape)

#     gained_df2,lost_df2,const_df2 = CTCF_TALL_modules.return_specific_binding_df(cancertype)
#     const_df = const_df.loc[const_df2.index]#;print(const_df.shape);exit()
    return gained_df,lost_df,const_df.iloc[:1000,:]

        

def main_cancertype(cancertype):
    
    indir = 'f3_compr_ttest_csv'
    outdir = 'f4_scatter_with_boxplot_figs'
    os.makedirs(outdir,exist_ok=True)
    viewregions = [20000,50000,100000,200000,500000,1000000,2000000]
    viewregions = [500000]
    
    for normalization in  ['raw']:
        for resolution in [5000]:
            for viewregion in viewregions:
                if viewregion > resolution and viewregion%resolution==0:
                    gained_df,lost_df,const_df = read_ttest_df(indir,cancertype,normalization,resolution,viewregion)
                    print(cancertype,normalization,resolution,viewregion,gained_df.shape,lost_df.shape,const_df.shape)
                    #print(gained_df);exit()
                    scatter_plot(gained_df,lost_df,const_df,cancertype,normalization,resolution,viewregion,outdir)
        
                    gained_total,gained_pos,gained_minus,gained_values= return_sig_nums(cancertype,gained_df);#print(gained_total,gained_pos,gained_minus)
                    lost_total,lost_pos,lost_minus,lost_values = return_sig_nums(cancertype,lost_df);#print(lost_total,lost_pos,lost_minus)
                    const_total,const_pos,const_minus,const_values = return_sig_nums(cancertype,const_df) ;#print(const_total,const_pos,const_minus)        

                    figname = outdir+os.sep+'{}_{}_res{}_view{}_box_bar.pdf'.format(cancertype,normalization,resolution,viewregion)
                    box_plot(cancertype,gained_values,lost_values,const_values,normalization,resolution,viewregion,figname)
                    #exit()

                
def main():

    cancertypes=['T-ALL','T-ALL-P1','T-ALL-P2','CRC']#,'PRAD','PRAD_TissueAdded']
    for cancertype in cancertypes:
        main_cancertype(cancertype)
        
#         exit()
               


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--viewregion', action = 'store', type = int,dest = 'viewregion', help = 'input file of', metavar = '<int>')
    parser.add_argument('-n', '--normalization', action = 'store', type = str,dest = 'normalization', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-c', '--chrom', action = 'store', type = str,dest = 'chrom', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
