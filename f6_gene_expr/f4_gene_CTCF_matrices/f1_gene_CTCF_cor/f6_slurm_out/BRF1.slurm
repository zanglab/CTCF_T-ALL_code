#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BRF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BRF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRI3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRI3BP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRICD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRINP1 
