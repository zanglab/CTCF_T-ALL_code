#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STX4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STX6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STX7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STX8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STXBP1 
