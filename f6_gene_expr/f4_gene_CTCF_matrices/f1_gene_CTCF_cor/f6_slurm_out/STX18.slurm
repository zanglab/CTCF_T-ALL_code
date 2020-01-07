#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STX18.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STX18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STX19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STX1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STX1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STX3 
