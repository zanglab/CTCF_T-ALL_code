#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCNF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNH1 
