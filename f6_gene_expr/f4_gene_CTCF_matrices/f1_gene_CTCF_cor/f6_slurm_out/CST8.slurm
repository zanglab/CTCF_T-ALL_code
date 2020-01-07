#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CST8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CST8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CST9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CST9L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSTA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSTB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSTF1 
