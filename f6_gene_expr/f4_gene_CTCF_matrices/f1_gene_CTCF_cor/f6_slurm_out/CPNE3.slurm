#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CPNE3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CPNE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPNE4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPNE5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPNE6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPNE7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPNE8 
