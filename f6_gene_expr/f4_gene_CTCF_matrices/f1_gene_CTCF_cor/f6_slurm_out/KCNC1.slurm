#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCNC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCND2 
