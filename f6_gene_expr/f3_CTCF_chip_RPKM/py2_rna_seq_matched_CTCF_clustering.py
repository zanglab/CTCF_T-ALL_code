import os,sys,argparse,glob,re,bisect
import numpy as np
import pandas as pd
from collections import Counter
import operator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#matplotlib.rcParams['agg.path.chunksize'] = 10000
matplotlib.rcParams['font.size']=12
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
#sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
from scipy.cluster.hierarchy import linkage,dendrogram,cut_tree


def hi_plot(df,k,labels,metric,figname):
    plt.figure(figsize=(5,25))
    Z = linkage(df,method = 'average',metric = metric)
    color_cut=Z[-(k - 1),2]
    dn = dendrogram(Z,labels = labels,color_threshold=color_cut,orientation='left',leaf_font_size=7)
    plt.savefig('{}_{}.pdf'.format(figname,metric),bbox_inches = 'tight',pad_inches = .1)
    plt.close()

             
def count_occurrence(num_thre):
    peak_occur_in_union_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/f4_lost_on_consistent_binding/f1_consistent_bindings/f0_infiles/union_summits_width_150_EachDataOverlapInfo.csv'
    #peak_occur_in_union_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/f4_lost_on_consistent_binding/f1_consistent_bindings/f0_infiles/union_summits_width_150_EachDataOverlapInfo_head10000.csv'
    peak_occur_in_union_df = pd.read_csv(peak_occur_in_union_file,sep="\t",index_col = 0) #;print(peak_occur_in_union_df);exit(0)
    all_sum = peak_occur_in_union_df.sum(axis=1);#print(all_sum);
    return all_sum[all_sum>num_thre].index #num_thre = 5, length=184688
    #return all_sum[(all_sum<584*(1-num_thre))&(all_sum>584*num_thre)].index #.1 ->66692 # .05 ->(95333, 584)

def rna_seq_matched_index():
    rna_seq_matched_file = 'geneExpr_matched_CTCF.xlsx'
    rna_seq_matched_df = pd.read_excel(rna_seq_matched_file,sep="\t",index_col = 0) #;print(rna_seq_matched_df.index);exit(0)
    return rna_seq_matched_df.index 
 
def get_signal_df():
    signal_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/f3_union_CTCF_regions/f3_signal_on_union_regions/f2_signals_on_union_bindings/signals_RPKM_on_all_CTCF_bindings_QuantileNormalized.csv'
    #signal_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/f3_union_CTCF_regions/f3_signal_on_union_regions/f2_signals_on_union_bindings/signals_RPKM_on_all_CTCF_bindings_QuantileNormalized_head10000.csv'
    signal_df = pd.read_csv(signal_file,sep="\t",index_col = 0)
    return signal_df

def get_matched_name_df():
    matched_name_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/fz_gene_CTCF_mete_info/CTCF_binding_clustering/all_CTCF_datasets_IF5_peak2000_manuRevised_20180805.csv'
    matched_name_df = pd.read_csv(matched_name_file,index_col=0)
    matched_name_df = matched_name_df.fillna('NA')
    #print(matched_name_df)
    return matched_name_df
    
def main():
    
    outdir = 'f2_rna_seq_matched_clustering_fig'
    os.makedirs(outdir,exist_ok=True)
    
    signal_df = get_signal_df();print(signal_df.shape)
    filtered_index = count_occurrence(5)
    signal_df = signal_df.loc[filtered_index];print(signal_df.shape)
    signal_df = signal_df.dropna()
    
    #signal_sum = signal_df.sum()
    #print(signal_sum[signal_sum==0]) # Series([], dtype: float64) # in case some vector=[0], not valid for correlation
    print(signal_df.shape);# exit(0) 
   
    # assign new column names
    matched_name_df = get_matched_name_df()
    filtered_index = rna_seq_matched_index()
    ori_colums = []
    new_column = []
    for column in signal_df.columns:
        if re.search('GSM',column):
            scolumn = column.split('_')[0]
        else:
            scolumn = '_'.join(column.split('_')[:-2])
        assert scolumn in matched_name_df.index
        if scolumn in filtered_index:
            #print(scolumn)
            ori_colums.append(column)
            a = matched_name_df.loc[scolumn,'cellline']
            b = matched_name_df.loc[scolumn,'tissue']
            c = int(matched_name_df.loc[scolumn,'peak_nums'])
            new_column.append('_'.join([scolumn,a,b,'Peak#'+str(c)]))
    
    signal_df = signal_df[ori_colums];#print(signal_df)
    signal_df.columns = new_column;#print(signal_df)
    
    #sns_clustering(signal_df,'{}/{}.pdf'.format(outdir,'sns_heatmap'))
    hi_plot(np.transpose(signal_df),10,signal_df.columns,'euclidean','{}/{}'.format(outdir,'hc'))
    hi_plot(np.transpose(signal_df),10,signal_df.columns,'correlation','{}/{}'.format(outdir,'hc'))




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
  
    main()
