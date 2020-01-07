#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AEN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AEN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AES 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AFAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AFAP1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AFAP1L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AFF1 
