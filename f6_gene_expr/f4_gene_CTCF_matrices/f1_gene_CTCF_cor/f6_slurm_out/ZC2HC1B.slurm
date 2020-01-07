#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZC2HC1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC2HC1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC2HC1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H11A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H12A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3H12B 
