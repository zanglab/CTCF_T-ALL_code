#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHRM4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRNA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRNA10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRNA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRNA3 
