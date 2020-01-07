#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NCOA5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NCOA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCOA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCOA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCOR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCOR1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCOR2 
