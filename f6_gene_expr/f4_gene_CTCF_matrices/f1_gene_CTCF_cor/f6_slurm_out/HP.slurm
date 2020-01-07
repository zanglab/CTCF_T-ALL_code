#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HP1BP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HPCA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HPCAL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HPCAL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HPD 
