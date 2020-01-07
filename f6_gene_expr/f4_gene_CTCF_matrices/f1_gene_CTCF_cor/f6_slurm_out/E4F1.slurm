#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o E4F1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i E4F1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EAF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EAPP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EARS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EBAG9 
