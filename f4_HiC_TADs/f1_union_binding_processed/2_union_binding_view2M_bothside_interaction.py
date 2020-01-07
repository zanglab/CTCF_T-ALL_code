
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
import CTCF_TALL_modules_new


chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']

all_bindings = CTCF_TALL_modules_new.return_occupancy_filtered()


##########################################
##########################################


def mirror_concat(df):

    for column_pos in np.arange(1,len(df.columns)):
        column = df.columns[column_pos]
        mir_column = '-{}'.format(column)
        new_values = np.append([0]*column_pos,df[column].values)[:len(df.index)]
        df[mir_column]=new_values
        #print(df[column],new_values);exit()
    return df
    

def write_out_binding_interactions_all_chroms(binding_df,data,normalization,viewregion,resolution,flag,outdir):
    
    binding_data_out = open(outdir+os.sep+'{}_{}_res{}_{}.csv'.format(data,normalization,resolution,flag),'w')
    
    columns = np.arange(-1*viewregion+resolution,viewregion,resolution)
    binding_data_out.write('{}\t{}\n'.format('id','\t'.join(map(str,columns))))
            
    for chrom in chroms:
        #chrom='chr22'
        print(chrom)
        matrix_file = 'f1_viewpoint_2M_interaction_abstraction/{}_{}_{}_res_{}_view_region_2000000.csv'.format(data,normalization,chrom,resolution)
        if os.path.isfile(matrix_file):
            with open(matrix_file) as matrix_inf:
                matrix_df = pd.read_csv(matrix_inf,sep='\t',index_col=0)
        #####
        # process the matrix df
        #####
        
            # normalization of the data, each interaction score is divided by average interaction score with same distances
            matrix_df = matrix_df/matrix_df.mean()#;print(matrix_file)
            matrix_df = mirror_concat(matrix_df)#;print(binding_data_df.loc[248900000][map(str,columns)].values);exit()
            binding_df_chr = binding_df.loc[binding_df['chr']==chrom]
            for binding_id in binding_df_chr.index:
                binding_pos = all_bindings.loc[binding_id].mid_position
                view_pos = binding_pos//resolution*resolution
                view_pos_binding_data = matrix_df.loc[view_pos][map(str,columns)].values  
                view_pos_binding_data = np.round(view_pos_binding_data,2)
                #print('\t'.join(map(str,view_pos_binding_data)))
                #print(view_pos_binding_data,len(columns),len(view_pos_binding_data));exit()
                binding_data_out.write('{}\t{}\n'.format(binding_id,'\t'.join(map(str,view_pos_binding_data))))
                #print(view_pos_binding_data);exit()
    binding_data_out.close() #;exit() 
     
     
     
    
        
        
def main(view_region,normalization,chrom):
    
    outdir = 'f2_union_binding_view2M_bothside_interaction'
    os.makedirs(outdir,exist_ok=True)
    
    viewregion = 2000000
    data_types =  ['Jurkat','A6010','Cutll1','PD31','PD9','HCT116','trans_colon1','trans_colon2','LNCaP','PrEC']
    
    for normalization in ['raw']:
        for data in data_types:
            for resolution in [5000]:          
                write_out_binding_interactions_all_chroms(all_bindings,data,normalization,viewregion,resolution,'union',outdir)







if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--viewregion', action = 'store', type = int,dest = 'viewregion', help = 'input file of', metavar = '<int>')
    parser.add_argument('-n', '--normalization', action = 'store', type = str,dest = 'normalization', help = 'input file of', metavar = '<str>')
    parser.add_argument('-c', '--chrom', action = 'store', type = str,dest = 'chrom', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.viewregion,args.normalization,args.chrom)
