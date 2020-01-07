#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SAXO2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SAXO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAYSD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SBDS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SBDSP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SBF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SBF1P1 
