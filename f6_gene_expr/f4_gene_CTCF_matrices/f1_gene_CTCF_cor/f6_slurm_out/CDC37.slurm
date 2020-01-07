#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDC37.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC37L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC42BPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC42BPB 
