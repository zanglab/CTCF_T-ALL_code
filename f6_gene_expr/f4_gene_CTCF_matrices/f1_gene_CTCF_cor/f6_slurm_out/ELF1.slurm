#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ELF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ELF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELFN1 
