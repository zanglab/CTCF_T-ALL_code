#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MIF4GD.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MIF4GD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MILR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MINK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MINOS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MINOS1-NBL1 
