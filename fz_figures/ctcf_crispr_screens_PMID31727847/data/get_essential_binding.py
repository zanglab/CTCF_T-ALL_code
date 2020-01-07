import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
#sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
#sns.despine(offset=0, trim=True)




def main():

    outdir = 'bed'
    os.makedirs(outdir,exist_ok=True)

    # ==== read table of result and ctcf binding info
    target_file='supp_tables/pnas.1908155116.sd04.xlsx'
    target_df = pd.read_excel(target_file,index_col=3,header=None)
    target_df.columns = ['chrom','start','end']
    # ==== remove those positive control
    result_file='supp_tables/pnas.1908155116.sd05.xlsx'
    result_df = pd.read_excel(result_file,index_col=0)
    result_df = result_df[result_df['Type']=='Peak']
    
    for cancertype in ['LNCaP','T47D']:
        beta_col = '{}.beta'.format(cancertype)
        fdr_col = '{}.fdr'.format(cancertype)
        fdr_thre = 0.25
        # ==== keep those CTCF binding with the 10% lowest beta-score
        beta_df = result_df.sort_values(by=[beta_col],ascending=True).iloc[:int(result_df.shape[0]*.2)]
        beta_df = beta_df.join(target_df)
        beta_df = beta_df.dropna(how='any',axis=0)
        beta_df.insert(0,'id',beta_df.index)
        # ==== or keep CTCF binding with FDR<0.25
        fdr_df = result_df[(result_df[fdr_col]<fdr_thre)]
        fdr_df = fdr_df.join(target_df)
        fdr_df = fdr_df.dropna(how='any',axis=0)
        fdr_df.insert(0,'id',fdr_df.index)
        
        fdr_df['start'] = fdr_df['start'].astype(int)
        fdr_df['end'] = fdr_df['end'].astype(int)
        beta_df['start'] = beta_df['start'].astype(int)
        beta_df['end'] = beta_df['end'].astype(int)
        
        beta_df.to_csv(outdir+os.sep+'{}_beta_filtered.csv'.format(cancertype))
        fdr_df.to_csv(outdir+os.sep+'{}_fdr_filtered.csv'.format(cancertype))
        
        beta_df[['chrom','start','end','id']].to_csv(outdir+os.sep+'{}_beta_filtered.bed'.format(cancertype),index=False,sep='\t',header=None)
        fdr_df[['chrom','start','end','id']].to_csv(outdir+os.sep+'{}_fdr_filtered.bed'.format(cancertype),index=False,sep='\t',header=None)
        combined_df = pd.concat([beta_df,fdr_df]).drop_duplicates()
        combined_df[['chrom','start','end','id']].to_csv(outdir+os.sep+'{}_fdr_combi_beta_filtered.bed'.format(cancertype),index=False,sep='\t',header=None)
        




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
