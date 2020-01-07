#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KLF3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF8 
