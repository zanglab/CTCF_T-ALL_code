import os,sys,argparse,glob
import numpy as np
import pandas as pd
import find_overlap_keep_info_NOT_sep_strand_asimport
import re
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=14
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
import CTCF_TALL_modules_new

#>>> set(df['annotation'].values)
#{'Exon', 'Promoter', "5' UTR", 'Distal', "3' UTR", 'Intron'}

annotation_marks = [ 'Intron','Exon', 'Promoter', 'Distal']

    
def compr_plot_box(annotation,flag,outdir):
    print(annotation)#;exit()
    plt.figure(figsize=(2.5,2.3))
    cancertypes = sorted(annotation.keys())
    cancertypes = annotation.keys()
    colors = ['lightgreen','cornflowerblue','lightgrey','lightsalmon'][::-1]
    labels = ["Promoter","Exon","Intron","Distal"][::-1]
    position = 1
    for cancertype in cancertypes:
        a = annotation[cancertype][labels[0]]
        b = annotation[cancertype][labels[1]]
        c = annotation[cancertype][labels[2]]
        d = annotation[cancertype][labels[3]]
        #e = annotation[cancertype][labels[4]]
        #f = annotation[cancertype][labels[5]]
        total = a+b+c+d
        print(a,b,c,d)
        g0 = plt.bar(position,100*a/total,bottom=0,width = .68, lw=0,color = colors[0],label = labels[0])
        g1 = plt.bar(position,100*b/total,bottom=100*a/total,width = .68, lw=0,color = colors[1],label = labels[1])
        g2 = plt.bar(position,100*c/total,bottom=100*(a+b)/total,width = .68, lw=0,color = colors[2],label = labels[2])
        g3 = plt.bar(position,100*d/total,bottom=100*(a+b+c)/total,width = .68, lw=0,color = colors[3],label = labels[3])
        #g4 = plt.bar(position,100*e/total,bottom=100*(a+b+c+d)/total,width = .5,color = colors[4],label = labels[4])
        #g5 = plt.bar(position,100*f/total,bottom=100*(a+b+c+d+e)/total,width = .5,color = colors[5],label = labels[5])
        position+=1
        
    #positions = [1,2,3]
    #g = plt.boxplot(box_vals,positions=positions,widths = .5,patch_artist=True)
    if flag=='lost':
        plt.ylabel('% CTCF binding'.format(flag),fontsize=16)
    plt.ylim([0,100])
    plt.xlim([0.3,6.7])
    plt.title('{}'.format(flag.capitalize()),fontsize=17)
    if flag=='gained':
        plt.legend([g3,g2,g1,g0],["Promoter","Exon","Intron","Intergenic"],loc=1,bbox_to_anchor=(1.7,1.05),fontsize=14,borderaxespad=0.,labelspacing=.2,handletextpad=0.2,handlelength=1,frameon=False)
    plt.axes().set_xticks([1.2,2.2,3.2,4.2,5.2,6.2]) 
    sns.despine(offset=0, trim=True)
    plt.axes().spines['bottom'].set_visible(False)
    ticklabels=['T-ALL','AML','BRCA','CRC','LUAD','PRAD']
    plt.axes().set_xticklabels(ticklabels,rotation=45,ha='right',fontsize=15) 
    plt.axes().tick_params(axis='x',direction='out', length=0, width=.8, colors='black')
    plt.savefig('{}/cancer_specific_{}_position_compr.pdf'.format(outdir,flag),bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)


def return_file_annotation(df):
    # for each file, return a diction of the annotation
    #print(df);exit()
    df_annotation_dic = {}
    df_annotation = df['annotation'].values
    for annotation in annotation_marks:
        df_annotation_dic[annotation] = len([i for i in df_annotation if i==annotation])
    
    return df_annotation_dic           


def main():
    
    outdir = 'f1_binding_location'
    os.makedirs(outdir,exist_ok=True)

    cancertypes=['T-ALL','AML','BRCA','CRC','LUAD','PRAD_TissueAdded']
    ticklabels=['T-ALL','AML','BRCA','CRC','LUAD','PRAD']
    
    gained_annotation = {}
    lost_annotation = {}
    const_annotation = {}
    for cancertype in cancertypes:
        gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype)
        gained_annotation[cancertype] = return_file_annotation(gained_df)
        lost_annotation[cancertype] = return_file_annotation(lost_df)
        
        compr_plot_box(gained_annotation,'gained',outdir)
        compr_plot_box(lost_annotation,'lost',outdir)





           
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--infile1', action = 'store', type = str,dest = 'infile1', help = 'input file to be compared/separated', metavar = '<file>')
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
