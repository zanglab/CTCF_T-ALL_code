#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EXO5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EXO5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOC3L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOC3L2 
