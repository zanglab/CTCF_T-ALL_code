#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HTR6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HTR6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTR7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTR7P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTRA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTRA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTRA3 
