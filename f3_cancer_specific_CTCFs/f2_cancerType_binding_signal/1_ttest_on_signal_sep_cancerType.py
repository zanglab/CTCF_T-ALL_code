import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
import scipy
from scipy import stats
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['agg.path.chunksize'] = 10000
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=1.1)
#sns.set_style("whitegrid", {'axes.grid' : False})

def fdr_adj_p(pvalues,p_index):
    # p_index here is only use to find the right column names
    df = pvalues.to_frame();#print(df);exit(0)
    n,k = len(pvalues),len(pvalues)      
    minimum = 1    
    for index in df.sort_values(by=p_index,ascending=False).index:
        pvalue = df.loc[index,p_index]
        fdr = n*pvalue/k  
        minimum = min(minimum,fdr)
        df.loc[index,'adj.p'] = minimum
        k=k-1     
    return df['adj.p']



def main(cancerType):

    outdir = 'f1_cancertype_signal_ttest_fdr'
    os.makedirs(outdir,exist_ok=True)
    signal_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f3_signals_on_union_CTCFs/f2_signals_on_union_bindings/signals_RPKM_on_all_CTCF_bindings_QuantileNormalized.csv'
#     signal_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f3_signals_on_union_CTCFs/f2_signals_on_union_bindings/signals_RPKM_on_all_CTCF_bindings_QuantileNormalized_head10000.csv'
    signals = pd.read_csv(signal_file,sep='\t',index_col=0)

    # read cancer and normal datasets/GSMID for each cancer type
    cancerType_data = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f0_infile/CancerTypes_datasets_201907.xlsx'
    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']

#     for cancerType in cancertypes:
    if 1:
        cancertype_df = pd.read_excel(cancerType_data,index_col=0,sheetname=cancerType)
        # cancer and normal GSM IDs
        cancercols = cancertype_df['CancerData'].dropna().values
        normalcols = cancertype_df['NormalData'].dropna().values
        print('\n',cancerType,'\tCancer-len:',len(cancercols),'\tNormal-len:',len(normalcols))
        print('\n',cancercols,'\n',normalcols,'\n')
        controlcols = [i for i in signals.columns if i not in cancercols]
        controlcols = [i for i in signals.columns if i not in normalcols]

        # incase you have to run the script again
#         if not os.path.isfile('{}/{}_signals_ttest_fdr.csv'.format(outdir,cancerType)):
        if 1:
            i=0
            ttest_file = open('{}/{}_signals_ttest.csv'.format(outdir,cancerType),'w')
            ttest_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('ID',\
            'cancer_mean','cancer_median','normal_mean','normal_median','ctrl_mean','ctrl_median',\
            'cancer_vs_other_stats','cancer_vs_other_pvalue','normal_vs_other_stats','normal_vs_other_pvalue',\
            'cancer_vs_normal_stats','cancer_vs_normal_pvalue'))
            for bindingID in signals.index:            
                binding_signal = signals.loc[bindingID]   
                cancer_sig = binding_signal[cancercols]
                normal_sig = binding_signal[normalcols]
                control_sig = binding_signal[controlcols]
                #print(cancer_sig.mean(),control_sig.mean());exit()
                cancer_vs_other_stats,cancer_vs_other_pvalue = stats.ttest_ind(cancer_sig,control_sig)
                normal_vs_other_stats,normal_vs_other_pvalue = stats.ttest_ind(normal_sig,control_sig)
                cancer_vs_normal_stats,cancer_vs_normal_pvalue = stats.ttest_ind(cancer_sig,normal_sig)
                #cancertype_pvalue.loc[bindingID,'stat_score'] = stats_value
                #cancertype_pvalue.loc[bindingID,'pvalue'] = pvalue
                ttest_file.write('{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.3e}\t{:.2f}\t{:.3e}\t{:.2f}\t{:.3e}\n'.format(bindingID,\
                cancer_sig.mean(),cancer_sig.median(),normal_sig.mean(),normal_sig.median(),control_sig.mean(),control_sig.median(),\
                cancer_vs_other_stats,cancer_vs_other_pvalue,normal_vs_other_stats,normal_vs_other_pvalue,cancer_vs_normal_stats,cancer_vs_normal_pvalue))
                #if i%10000==0:
                #    print(i,bindingID);#print(cancer_sig,control_sig)
                #i=i+1
            ttest_file.close()    
            cancertype_pvalue = pd.read_csv('{}/{}_signals_ttest.csv'.format(outdir,cancerType),sep='\t',index_col=0)
            cancertype_pvalue['cancer_vs_other_adjp'] =  fdr_adj_p(cancertype_pvalue['cancer_vs_other_pvalue'],'cancer_vs_other_pvalue')
            cancertype_pvalue['normal_vs_other_adjp'] =  fdr_adj_p(cancertype_pvalue['normal_vs_other_pvalue'],'normal_vs_other_pvalue')
            cancertype_pvalue['cancer_vs_normal_adjp'] =  fdr_adj_p(cancertype_pvalue['cancer_vs_normal_pvalue'],'cancer_vs_normal_pvalue')
            cancertype_pvalue.to_csv('{}/{}_signals_ttest_fdr.csv'.format(outdir,cancerType),sep='\t')
            #print(cancertype_pvalue);exit(0)
            


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of bed fromat, union all the overlapping regions', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
   # parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile)
