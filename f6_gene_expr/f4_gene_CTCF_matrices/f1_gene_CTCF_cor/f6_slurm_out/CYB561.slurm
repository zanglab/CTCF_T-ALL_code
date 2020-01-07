#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYB561.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYB561 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYB561A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYB561D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYB561D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYB5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYB5B 
