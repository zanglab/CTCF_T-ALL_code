#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNHIT3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNHIT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNHIT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNRD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNRF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNRF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNRF2P1 
