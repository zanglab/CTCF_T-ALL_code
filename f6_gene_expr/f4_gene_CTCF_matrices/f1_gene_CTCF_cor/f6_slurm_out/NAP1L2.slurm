#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NAP1L2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NAP1L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAP1L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAP1L4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAP1L5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAP1L6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAPA 
