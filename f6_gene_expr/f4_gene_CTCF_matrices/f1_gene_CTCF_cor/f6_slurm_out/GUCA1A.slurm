#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GUCA1A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GUCA1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUCA1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUCA1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUCA2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUCA2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUCD1 
