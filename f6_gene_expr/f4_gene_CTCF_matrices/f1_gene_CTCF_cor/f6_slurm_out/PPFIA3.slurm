#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPFIA3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPFIA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPFIA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPFIBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPFIBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPHLN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIA 
