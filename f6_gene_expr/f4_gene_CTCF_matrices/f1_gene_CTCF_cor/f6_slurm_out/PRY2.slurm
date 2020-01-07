#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRY2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRY2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSAPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSCA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSD 
