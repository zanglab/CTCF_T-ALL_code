#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GEMIN6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GEMIN6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GEMIN7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GEMIN8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GEMIN8P4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GEN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GET4 
