#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GSDMD.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GSDMD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSG1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSG1L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSK3A 
