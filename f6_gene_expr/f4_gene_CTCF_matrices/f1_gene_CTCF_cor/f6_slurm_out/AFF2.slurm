#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AFF2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AFF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AFF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AFF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AFG3L1P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AFG3L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AFM 
