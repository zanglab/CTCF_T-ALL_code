import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']

def write_out_juicer_format_matrix(order_index_file,matrix_file,outdir,data,resolution,flag):

    #order_index_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/11_HiC_analysis/f1_preprocess/HiC_Pro/HiC_fastq_split/A6010/A6010_out/hic_results/matrix/fastq/raw/20000/fastq_20000_abs.bed'
    #matrix_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/11_HiC_analysis/f1_preprocess/HiC_Pro/HiC_fastq_split/A6010/A6010_out/hic_results/matrix/fastq/raw/20000/fastq_20000.matrix'
    #matrix_file = 'test_matrix.txt'
    with open(order_index_file) as order_index_inf: #, open(matrix_file) as matrix_inf:
        index_df = pd.read_csv(order_index_inf,sep='\t',header=None)
        #matrix_df = pd.read_csv(matrix_inf,sep='\t',header=None)
    index_df.columns = ['chr','x','y','id']
    #matrix_df.columns = ['id_a','id_b','score']
    
    for chrom in chroms:
        outfile = open(outdir+os.sep+'{}_{}_{}_{}.matrix'.format(data,resolution,flag,chrom),'w')
        index_df_chr = index_df[index_df['chr']==chrom]
        index_df_chr.index = index_df_chr['id']
        index_dict = index_df_chr[['x']].to_dict()['x']
        #print(index_df_chr);exit()
        keys = index_dict.keys()
        
        matrix_inf = open(matrix_file)
        line = matrix_inf.readline()
        while line:
            sline = line.strip().split()
            id_a,id_b,score = int(sline[0]),int(sline[1]),float(sline[2])
            if (id_a in keys) and (id_b in keys):
                outfile.write('{}\t{}\t{:.2f}\n'.format(index_dict[id_a],index_dict[id_b],score))
            #print(sline);exit()
            line = matrix_inf.readline()
        matrix_inf.close()
        outfile.close()
        #print(index_dict);exit()

def main(data_type):

    par_outdir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f4_HiC_TADs/f1_union_binding_processed/f0_transformed_matrix'
    par_indir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f4_HiC_TADs/f0_HiC_process/prostate_HiC_Pro_BglII/prostate_bglii/hic_results/matrix'
    
    #for data_type in ['A6010','Cutll1','Jurkat','PD31','PD9']:
    if data_type in ['LNCaP','PrEC']:
        for resolution in [5000]:
            if 1:
                outdir = par_outdir+os.sep+'{}/{}'.format(data_type,resolution)
                os.makedirs(outdir,exist_ok=True)
                index_file = par_indir+os.sep+'{}/raw/{}/{}_{}_abs.bed'.format(data_type,resolution,data_type,resolution)
                raw_matrix = par_indir+os.sep+'{}/raw/{}/{}_{}.matrix'.format(data_type,resolution,data_type,resolution)
#                 iced_matrix = par_indir+os.sep+'{}/iced/{}/{}_{}_iced.matrix'.format(data_type,resolution,data_type,resolution)
                print(data_type,resolution,os.path.isfile(index_file),os.path.isfile(raw_matrix))
                write_out_juicer_format_matrix(index_file,raw_matrix,outdir,data_type,resolution,'raw')
#                 write_out_juicer_format_matrix(index_file,iced_matrix,outdir,data_type,resolution,'iced')
    
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--data', action = 'store', type = str,dest = 'data', help = 'input file of', metavar = '<int>')
    #parser.add_argument('-r', '--resolution', action = 'store', type = int,dest = 'resolution', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.data)
