import os,sys,argparse,glob
import numpy as np
import pandas as pd
import find_overlap_keep_info_NOT_sep_strand_asimport


log_sig_file='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/15_TCGA_patient_data/ATAC_seq/data/Count_matrices/TCGA-ATAC_PanCan_Log2Norm_Counts.txt'
# log_sig_file='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/15_TCGA_patient_data/ATAC_seq/data/Count_matrices/TCGA-ATAC_PanCan_Log2Norm_Counts_head10000.txt'
raw_sig_file='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/15_TCGA_patient_data/ATAC_seq/data/Count_matrices/TCGA-ATAC_PanCan_Raw_Counts.txt'
# raw_sig_file='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/15_TCGA_patient_data/ATAC_seq/data/Count_matrices/TCGA-ATAC_PanCan_Raw_Counts_head10000.txt'
mynorm_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/15_TCGA_patient_data/ATAC_seq/data/Count_matrices/mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.txt"


def return_binding_signal(sig_file,flag,outdir):

    outdir = outdir+os.sep+flag
    os.makedirs(outdir,exist_ok=True)
    
    #with open(sig_file) as sig_inf:
    #    sig_df = pd.read_csv(sig_inf,sep='\t')
    
    pyfile="find_overlap_keep_info_NOT_sep_strand_revised.py"
    bed_file_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/fu_feature_combination/f3_cancer_specific_binding"
    binding_files = glob.glob(bed_file_dir+os.sep+'*.bed')
    for binding_file in binding_files:
        basename = os.path.basename(binding_file).split('.bed')[0];print('\n## ==== ',basename)
        outfile = '{}/{}_ATAC_sig.bed'.format(outdir,basename)
        print('python {} -a {} -b {} -s hg38 -p {} -e2 0'.format(pyfile,sig_file,binding_file,outfile))
#         os.system('python {} -a {} -b {} -s hg38 -p {} -e2 0'.format(pyfile,sig_file,binding_file,outfile))
#         exit()


#### ==== constitutive bindings ====
    constitutive_file='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_constitutive_0.8.bed'
    outfile = '{}/constitutive_ATAC_sig.bed'.format(outdir)
    print('\n## === constitutive')
    print('python {} -a {} -b {} -s hg38 -p {} -e2 0'.format(pyfile,sig_file,constitutive_file,outfile))
#     os.system('python {} -a {} -b {} -s hg38 -p {} -e2 0'.format(pyfile,sig_file,binding_file,outfile))



def main():

    outdir = 'f1_ATAC_sig_abstract'
    os.makedirs(outdir,exist_ok=True)
    
    #return_binding_signal(raw_sig_file,'raw',outdir)
    #return_binding_signal(log_sig_file,'log2norm',outdir)
    return_binding_signal(mynorm_file,'mynorm',outdir)
    
    #a = pd.read_csv('f1_abstract_sig/raw/BLCA_basicCTCF_overlapped_Breast_gained_ATAC_sig.bed',sep='\t')
    #print(a.iloc[3,8])
           
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--infile1', action = 'store', type = str,dest = 'infile1', help = 'input file to be compared/separated', metavar = '<file>')
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
