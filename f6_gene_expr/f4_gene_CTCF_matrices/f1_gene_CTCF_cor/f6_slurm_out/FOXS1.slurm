#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FOXS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FPGS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FPGT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FPGT-TNNI3K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FPR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FPR2 
