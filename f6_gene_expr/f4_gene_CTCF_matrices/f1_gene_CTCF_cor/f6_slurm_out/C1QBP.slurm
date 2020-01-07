#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C1QBP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QL4 
