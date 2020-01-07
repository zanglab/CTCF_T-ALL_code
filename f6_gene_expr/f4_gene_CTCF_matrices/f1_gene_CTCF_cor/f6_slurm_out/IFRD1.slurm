#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IFRD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IFRD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFRD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFT122 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFT140 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFT172 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFT20 
