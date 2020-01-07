#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CKAP5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CKAP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CKB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CKLF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CKLF-CMTM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CKM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CKMT1A 
