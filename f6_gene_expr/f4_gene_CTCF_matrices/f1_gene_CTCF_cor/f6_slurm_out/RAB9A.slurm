#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAB9A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB9A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB9B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB9BP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RABAC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RABEP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RABEP2 
