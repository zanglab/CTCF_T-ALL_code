#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o THBS4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i THBS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THEG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THEGL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THEM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THEM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THEM6 
