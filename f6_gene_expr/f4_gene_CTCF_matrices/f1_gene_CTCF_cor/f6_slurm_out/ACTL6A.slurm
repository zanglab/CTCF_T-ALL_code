#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACTL6A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTL6A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTL6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTL7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTL7B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTL8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTL9 
