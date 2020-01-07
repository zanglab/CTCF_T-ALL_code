#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMIM19.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM23 
