#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCND3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCND3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNE4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNE5 
