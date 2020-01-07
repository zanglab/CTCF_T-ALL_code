#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZSCAN22.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN30 
