#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCNJ6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNK10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNK12 
