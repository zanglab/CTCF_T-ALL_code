#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAD2L2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAD2L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MADCAM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MADD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAEA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAEL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAF 
