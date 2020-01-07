import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=14
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"



chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']

cancer_type_hic_filename = {'T-ALL':['Jurkat','A6010'],\
                                'T-ALL-P1':['PD9','A6010'],\
                                'T-ALL-P2':['PD31','A6010'],\
                                'CRC':['HCT116','trans_colon1'],\
                                'PRAD':['LNCaP','PrEC'],\
                                'PRAD_TissueAdded':['LNCaP','PrEC']}

##########################################
##########################################


def collect_interaction_compr_all_chroms(cancertype,data_type,normalization,resolution,viewregion):
        
    df_all = pd.DataFrame()
    cancer_file = 'f1_gained_lost_binding_interactions'+os.sep+'{}_{}_ctcf_interactions_in_{}_{}_res{}.csv'.format(cancertype,data_type,cancer_type_hic_filename[cancertype][0],normalization,resolution)
    normal_file = 'f1_gained_lost_binding_interactions'+os.sep+'{}_{}_ctcf_interactions_in_{}_{}_res{}.csv'.format(cancertype,data_type,cancer_type_hic_filename[cancertype][1],normalization,resolution)
    
    columns = np.arange(-1*viewregion+resolution,viewregion,resolution)
    with open(cancer_file) as cancer_inf, open(normal_file) as normal_inf:
        cancer_df = pd.read_csv(cancer_inf,index_col=0)
        normal_df = pd.read_csv(normal_inf,index_col=0)
    
    for index in cancer_df.index:
        cancer_sig = cancer_df.loc[index][map(str,columns)].values#;print(cancer_sig);exit()
        normal_sig = normal_df.loc[index][map(str,columns)].values
        stats_score,pvalue = stats.ttest_rel(cancer_sig,normal_sig)
        log_stats_score,log_pvalue = stats.ttest_rel(np.log2(cancer_sig+0.01),np.log2(normal_sig+0.01))

        df_all.loc[index,'stats'] = stats_score
        df_all.loc[index,'pvalue'] = pvalue
        df_all.loc[index,'log_stats'] = log_stats_score
        df_all.loc[index,'log_pvalue'] = log_pvalue
        df_all.loc[index,'{}_median'.format(cancer_type_hic_filename[cancertype][0])] = np.median(cancer_sig)
        df_all.loc[index,'{}_median'.format(cancer_type_hic_filename[cancertype][1])] = np.median(normal_sig)
        df_all.loc[index,'{}_mean'.format(cancer_type_hic_filename[cancertype][0])] = np.mean(cancer_sig)
        df_all.loc[index,'{}_mean'.format(cancer_type_hic_filename[cancertype][1])] = np.mean(normal_sig)
    #print(df_all);exit()
    return df_all





def main():

    outdir = 'f3_compr_ttest_csv'
    os.makedirs(outdir,exist_ok=True)

    cancertypes=['T-ALL','T-ALL-P1','T-ALL-P2','CRC','PRAD','PRAD_TissueAdded']
    for cancertype in cancertypes:
        for normalization in ['raw']:
            for resolution in [5000]:
                for viewregion in [20000,50000,100000,200000,500000,1000000,2000000]:
                    if viewregion>resolution and viewregion%resolution==0:
                        for data_type in ['gained','lost','ctrl']:
                            out_df = collect_interaction_compr_all_chroms(cancertype,data_type,normalization,resolution,viewregion)
                            out_df.to_csv(outdir+os.sep+'{}_{}_{}_res-{}_viewregion-{}_paired_ttest.csv'.format(cancertype,data_type,normalization,resolution,viewregion))
                            #exit()
            
        



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--resolution', action = 'store', type = int,dest = 'resolution', help = 'input file of', metavar = '<int>')
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
