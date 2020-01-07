#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o YIPF5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i YIPF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YIPF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YIPF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YJEFN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YKT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YLPM1 
