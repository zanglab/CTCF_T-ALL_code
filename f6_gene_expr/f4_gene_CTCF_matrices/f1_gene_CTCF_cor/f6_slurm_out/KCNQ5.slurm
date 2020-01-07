#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCNQ5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNQ5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNRG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNT1 
