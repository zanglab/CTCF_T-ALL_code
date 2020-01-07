#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZDHHC17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC21 
