#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MTMR3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MTMR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTMR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTMR6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTMR7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTMR8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTMR9 
