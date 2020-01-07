#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GLIS3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GLIS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLMN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLMP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLOD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLOD5 
