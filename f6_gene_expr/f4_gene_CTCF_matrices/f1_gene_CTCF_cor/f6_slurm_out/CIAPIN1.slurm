#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CIAPIN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CIAPIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CIART 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CIB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CIB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CIB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CIB4 
