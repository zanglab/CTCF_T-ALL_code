#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCNH8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNH8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNIP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNIP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNJ1 
