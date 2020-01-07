#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CNPY1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CNPY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNPY2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNPY3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNPY4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNR2 
