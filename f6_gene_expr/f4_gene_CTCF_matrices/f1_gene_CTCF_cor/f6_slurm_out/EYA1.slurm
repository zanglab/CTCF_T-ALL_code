#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EYA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EYA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EYA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EYA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EYA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EYS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EZH1 
