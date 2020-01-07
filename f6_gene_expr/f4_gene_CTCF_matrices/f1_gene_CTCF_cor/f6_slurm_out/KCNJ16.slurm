#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCNJ16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ5 
