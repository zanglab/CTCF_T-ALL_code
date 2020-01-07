#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TPCN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TPCN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPD52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPD52L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPD52L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPD52L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPGS1 
