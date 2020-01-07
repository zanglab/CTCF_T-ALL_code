#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATE1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATE1-AS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATF4 
