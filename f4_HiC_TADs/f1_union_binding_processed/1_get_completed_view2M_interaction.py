import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect

##########################################
## basic models for commonly used data 
##########################################
chrom_size_file = '/nv/vol190/zanglab/zw5j/data/Genome/hg38/hg38_clean.chrom.sizes'
chrom_size_df = pd.read_csv(chrom_size_file,sep='\t',header=None,index_col=0);
chrom_size_df.columns = ['len']

chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']

#union_CTCF_kept_ids_sample_thre_2_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/fz_gene_CTCF_mete_info/commonly_used_modules_files/union_CTCF_kept_ids_sample_thre_2.csv'
#all_binding_middle_domain_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/fz_gene_CTCF_mete_info/commonly_used_modules_files/all_CTCF_domainInfo.csv'
#with open(union_CTCF_kept_ids_sample_thre_2_file) as filter_inf, open(all_binding_middle_domain_file) as all_inf:
#    kept_ids = [int(i.strip()) for i in filter_inf.readlines()]
#    all_bindings = pd.read_csv(all_inf,sep='\t')

#gained_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/fz_combine_all_features/f1_combined_csv/T-ALL_CTCF_gained_AllFeatures.csv'
#lost_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/fz_combine_all_features/f1_combined_csv/T-ALL_CTCF_lost_AllFeatures.csv'
#constitutive_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/f4_lost_on_consistent_binding/f1_consistent_bindings/f1_consistent_bindings/constitutive_CTCF_bindings_thre0.8.bed'

#with open(gained_file) as gained_inf, open(lost_file) as lost_inf, open(constitutive_file) as constitutive_inf:
#    gained_df = pd.read_csv(gained_inf,index_col=3)
#    lost_df = pd.read_csv(lost_inf,index_col=3)
#    const_df = pd.read_csv(constitutive_inf,sep='\t',index_col = 3,header=None)
#    const_df.columns = ['chr','start','end','score','strand']

##########################################
##########################################



def write_out_interactions(matrix_df,view_region,outdir,data,resolution,normalization,chrom):
    
    # for each position/view_pos in chrom, get the interaction with downstream loci within view_region
    outfile_name = outdir+os.sep+'{}_{}_{}_res_{}_view_region_{}.csv'.format(data,normalization,chrom,resolution,view_region)
    if os.path.isfile(outfile_name):
        print('Do not re-write the files')
        exit()
    out_file = open(outfile_name,'w')
    view_poses = np.arange(0,chrom_size_df.loc[chrom,'len'],resolution);print(chrom,'len-view_poses:\t',len(view_poses))
    out_file.write('{}\t{}\n'.format('dis','\t'.join(map(str,np.arange(0,view_region,resolution)))))
    for view_pos in view_poses:
        #view_pos = 50785000
        view_df = matrix_df[matrix_df['x']==view_pos]
        view_df = view_df[['score','dis']]; 
        #print(view_pos,view_df)  ;exit()     
        standard_df = pd.DataFrame(columns = ['dis'],index = np.arange(0,view_region,resolution))
        standard_df['dis']=standard_df.index
        # the matrix_df has been restricted within view-region, here we could use outer or left
        # here is to get each of the interaction scores within the view region, if the interaction score exist in view_df, fill with 0
        merge_df = pd.merge(standard_df,view_df,how='outer').fillna(0)#;print(merge_df);exit()
        out_file.write('{}\t{}\n'.format(view_pos,'\t'.join(map(str,merge_df.score.values))));
    out_file.close()
    


def return_transformed_matrix_dir(data):
    if data in ['Jurkat','A6010','Cutll1','PD31','PD9']:
        indir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/11_HiC_analysis/f1_preprocess/union_binding_processed_data/panos/transformed_matrix'
    elif data in ['HCT116','trans_colon1','trans_colon2']:
        indir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/11_HiC_analysis/f1_preprocess/union_binding_processed_data/Public/transformed_matrix'
    elif data in ['LNCaP','PrEC']:
        indir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f4_HiC_TADs/f1_union_binding_processed/f0_transformed_matrix'   
    return indir

def main(chrom):
    
    outdir = 'f1_viewpoint_2M_interaction_abstraction'
    os.makedirs(outdir,exist_ok=True)
    view_region = 2000000

    for data in ['Jurkat','A6010','Cutll1','PD31','PD9','HCT116','trans_colon1','trans_colon2','LNCaP','PrEC']:
        indir = return_transformed_matrix_dir(data)
        for resolution in [5000]:
            for normalization in ['raw']:
                matrix_file = indir+os.sep+'{}/{}/{}_{}_{}_{}.matrix'.format(data,resolution,data,resolution,normalization,chrom)
                print(data,resolution,normalization,os.path.isfile(matrix_file))
                with open(matrix_file) as matrix_inf:
                    matrix_df = pd.read_csv(matrix_inf,sep='\t',header=None)
                    matrix_df.columns = ['x','y','score'.format(data)]
                matrix_df = matrix_df[matrix_df['y']-matrix_df['x']<view_region]
                # get the distances between x and y
                matrix_df['dis'] = matrix_df['y']-matrix_df['x']
                #print(matrix_df)
                write_out_interactions(matrix_df,view_region,outdir,data,resolution,normalization,chrom)
                #exit()
    




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-v', '--viewregion', action = 'store', type = int,dest = 'viewregion', help = 'input file of', metavar = '<int>')
    #parser.add_argument('-d', '--data', action = 'store', type = str,dest = 'data', help = 'input file of', metavar = '<int>')
    #parser.add_argument('-n', '--normalization', action = 'store', type = str,dest = 'normalization', help = 'input file of', metavar = '<str>')
    parser.add_argument('-c', '--chrom', action = 'store', type = str,dest = 'chrom', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.chrom)
