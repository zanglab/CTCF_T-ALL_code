#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IGDCC3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IGDCC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGDCC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGF1R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGF2BP1 
