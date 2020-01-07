#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LY6G5B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LY6G5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LY6G5C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LY6G6C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LY6G6D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LY6G6E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LY6G6F 
