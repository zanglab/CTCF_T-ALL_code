#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SGMS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SGMS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGPP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGPP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGSH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGSM1 
