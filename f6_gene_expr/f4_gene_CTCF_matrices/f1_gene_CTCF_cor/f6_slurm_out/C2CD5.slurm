#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C2CD5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C2CD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf27A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf42 
