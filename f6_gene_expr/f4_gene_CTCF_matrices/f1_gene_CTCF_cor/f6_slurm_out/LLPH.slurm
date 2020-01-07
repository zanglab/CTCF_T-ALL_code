#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LLPH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LLPH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMAN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMAN1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMAN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMAN2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMBR1 
