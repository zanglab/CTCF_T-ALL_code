#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KLC3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KLC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF12 
