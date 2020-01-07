#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IL4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL4I1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL4R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL5RA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL6 
