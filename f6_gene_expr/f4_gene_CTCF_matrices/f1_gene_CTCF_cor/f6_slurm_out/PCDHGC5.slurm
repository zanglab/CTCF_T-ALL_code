#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCDHGC5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCED1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCED1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCF11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCGF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCGF2 
