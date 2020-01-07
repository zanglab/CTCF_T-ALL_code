#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IFIT5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IFIT5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFITM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFITM10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFITM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFITM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFITM4P 
