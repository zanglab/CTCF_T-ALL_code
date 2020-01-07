#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPP14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPP21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPP25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPP25L 
