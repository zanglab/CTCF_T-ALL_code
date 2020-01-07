#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCNK4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNK5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNK6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNK7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNK9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNMA1 
