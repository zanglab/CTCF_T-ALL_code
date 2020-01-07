#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCNJ10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ15 
