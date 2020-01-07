#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRFAP1L1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRFAP1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRGBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRGPRD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRGPRE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRGPRF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRGPRG 
