#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POLR1A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR1D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR1E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2A 
