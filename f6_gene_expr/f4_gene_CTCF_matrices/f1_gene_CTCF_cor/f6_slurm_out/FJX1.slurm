#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FJX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FJX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP1A 
