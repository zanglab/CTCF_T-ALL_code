import os,sys,argparse,glob
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport


cancertypes = ['BRCA','CRC','LUAD','PRAD','PRAD_TissueAdded']
cacnertype_TCGA_matchname = {'BRCA':'BRCA','CRC':'COAD','LUAD':'LUAD','PRAD':'PRAD','PRAD_TissueAdded':'PRAD'}



def return_binding_signal(flag,outdir):

    outdir = outdir+os.sep+flag
    os.makedirs(outdir,exist_ok=True)
        
    pyfile="find_overlap_keep_info_NOT_sep_strand_revised.py"
    bed_file_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/fu_feature_combination/f3_cancer_specific_binding"
    norm_matrix_file_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f5_TCGA_ATAC/f0_cancerType_norm_matrix/overlapped"
    
    for cancertype in cancertypes:
        sig_file = norm_matrix_file_dir+os.sep+'ATAC_PanCan_Log2_QN_{}_overlapped.txt'.format(cacnertype_TCGA_matchname[cancertype])
        for ctcf_type in ['gained','lost']:
            binding_file = bed_file_dir+os.sep+'{}_{}.bed'.format(cancertype,ctcf_type)
            basename = os.path.basename(binding_file).split('.bed')[0];print('\n## ==== ',basename)
            outfile = '{}/{}_ATAC_sig.bed'.format(outdir,basename)
            print('python {} -a {} -b {} -s hg38 -p {} -e2 0'.format(pyfile,sig_file,binding_file,outfile))
#         os.system('python {} -a {} -b {} -s hg38 -p {} -e2 0'.format(pyfile,sig_file,binding_file,outfile))
#         exit()


#### ==== constitutive bindings ====
        constitutive_file='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_constitutive_0.8.bed'
        outfile = '{}/{}_constitutive_ATAC_sig.bed'.format(outdir,cancertype)
        print('\n## === constitutive')
        print('python {} -a {} -b {} -s hg38 -p {} -e2 0'.format(pyfile,sig_file,constitutive_file,outfile))
#     os.system('python {} -a {} -b {} -s hg38 -p {} -e2 0'.format(pyfile,sig_file,binding_file,outfile))



def main():

    outdir = 'f1_ATAC_sig_abstract'
    os.makedirs(outdir,exist_ok=True)
    
    return_binding_signal('mynorm',outdir)
    



           
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--infile1', action = 'store', type = str,dest = 'infile1', help = 'input file to be compared/separated', metavar = '<file>')
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
