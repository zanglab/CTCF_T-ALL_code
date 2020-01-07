#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANKHD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKHD1-EIF4EBP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKIB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKLE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKLE2 
