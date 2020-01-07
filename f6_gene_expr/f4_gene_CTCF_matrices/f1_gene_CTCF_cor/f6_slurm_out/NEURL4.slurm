#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NEURL4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NEURL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEUROD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEUROD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEUROD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEUROG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEUROG2 
