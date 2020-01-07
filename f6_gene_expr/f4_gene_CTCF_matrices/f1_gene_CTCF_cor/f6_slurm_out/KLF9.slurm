#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KLF9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHDC10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHDC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHDC4 
