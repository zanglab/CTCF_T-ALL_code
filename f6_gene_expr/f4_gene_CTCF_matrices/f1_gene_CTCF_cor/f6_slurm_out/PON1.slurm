#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PON1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PON1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PON2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PON3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POP5 
