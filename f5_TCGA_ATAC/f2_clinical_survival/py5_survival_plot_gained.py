import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=18
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.5)

sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style('ticks')
import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
# from stats_pvalues import irwin_hall_cdf
from lifelines.statistics import logrank_test
#from lifelines.plotting import add_at_risk_counts
from lifelines import KaplanMeierFitter


def survival_for_two(df,treat,ctrl,legends,title,figname):
    
    # select the time and status info for treat and control group
    ix = df['group'] == treat
    
    t1 = df.loc[ix]['time']
    e1 = df.loc[ix]['status'] 
    t2 = df.loc[~ix]['time']
    e2 = df.loc[~ix]['status']
    
    results = logrank_test(t1,t2,event_observed_A = e1,event_observed_B = e2)
    pvalue = results.p_value;print('pvalue:\t{}'.format(pvalue))
    
    # survival curves
    plt.figure(figsize=(3.,3.))
    ax = plt.subplot(111)
    
    kmf_control = KaplanMeierFitter()
    #g1 = kmf_control.fit(t1, e1, label=legends[0]).plot(ax=ax,show_censors=True,\
    g1 = kmf_control.fit(t1, e1).plot(ax=ax,show_censors=True,\
                        censor_styles={'ms': 12, 'marker': '+'},ci_show=False,c='red',ls='-')
    
    kmf_exp = KaplanMeierFitter()
    #g2 = kmf_exp.fit(t2, e2, label=legends[1]).plot(ax=ax,show_censors=True,\
    g2 = kmf_exp.fit(t2, e2).plot(ax=ax,show_censors=True,\
                    censor_styles={'ms': 12, 'marker': '+'},ci_show=False,c='k',ls='--')
    
    handles, labels = ax.get_legend_handles_labels();print(labels)
    lg = ax.legend(handles[1::2], legends,loc='lower left',borderaxespad=-.15,handletextpad=.2,labelspacing=.3,handlelength=1,frameon=False)
    if pvalue<0.2:
         plt.axes().text(df['time'].max()*0.45,0.45,'p={:.2f}'.format(pvalue),fontsize=16,ha='center')
    plt.ylim([-0.02,1.05])
#     plt.xlim([0,max_val*1])
    plt.title(title,fontsize=22)
    plt.xlabel('Days',fontsize=22)
    plt.ylabel('Survival probability',fontsize=22)
    plt.savefig(figname,bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
    plt.close()
    return results



def return_survival_df_50th(clinical_df):

    # get the clinical info, using top/bottom 25% as cutoff
    a_index = clinical_df.index[:int(len(clinical_df.index)*.5)]
    b_index = clinical_df.index[-1*int(len(clinical_df.index)*.5):]
    print('treat:\t',len(a_index),'\nctrl:\t',len(b_index))
    df = pd.DataFrame(index = a_index.union(b_index))
    df.loc[a_index,'group']='treat'   
    df.loc[b_index,'group']='ctrl'  
    df['status']= clinical_df.loc[df.index]['vital_status'] 
#     df.loc[df['status']=='Dead','status']=1
#     df.loc[df['status']=='Alive','status']=0;print(df.shape)
    df.loc[df['status']=='dead','status']=1
    df.loc[df['status']=='alive','status']=0;print(df.shape)
    df['time'] = clinical_df[['days_to_death','days_to_last_follow_up']].replace('--',0).astype(float).sum(axis=1)   
#     df['time'] = clinical_df[['days_to_death']].replace('--',0).astype(float).sum(axis=1)   
    return df

   
def main():

    outdir = 'f5_survival_figs_gained'
    os.makedirs(outdir,exist_ok=True)
    
    combined_file_dir='f3_combined_ATAC_sig_diff_score_and_clinical_info'
    cancertypes = ['BRCA','CRC','LUAD','PRAD','PRAD_TissueAdded']
    cancertypes = ['BRCA','CRC','LUAD']
    name_types_cancer_match={'BRCA':'BRCA','CRC':'COAD','LUAD':'LUAD','PRAD':'PRAD','PRAD_TissueAdded':'PRAD'}

    for name_type in cancertypes:
        combined_file = combined_file_dir+os.sep+'{}_gained_irwin_hall_with_clinical_info.csv'.format(name_type)
        clinical_df = pd.read_csv(combined_file,sep=',',index_col=0) 
        #print(clinical_df);exit()
        treat='treat'
        ctrl='ctrl'
        legends = ['more open','less open']
        title=name_types_cancer_match[name_type]
        print('\n\n=====\n',title)

        ### separate patient into two groups evenly, rank by mean different score
        print('cutoff: top/bottom 50%')
        df = return_survival_df_50th(clinical_df)
        df.to_csv(outdir+os.sep+'{}_survival_50th.csv'.format(name_type))
        figname=outdir+os.sep+'{}_survival_50th.pdf'.format(name_type)
        survival_for_two(df,treat,ctrl,legends,title,figname)




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
