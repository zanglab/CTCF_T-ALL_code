#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PIGS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGV 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGW 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGX 
