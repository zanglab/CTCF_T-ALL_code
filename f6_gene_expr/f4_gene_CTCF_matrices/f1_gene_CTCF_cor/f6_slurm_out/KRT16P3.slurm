#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRT16P3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT16P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT17P5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT18P55 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT19 
