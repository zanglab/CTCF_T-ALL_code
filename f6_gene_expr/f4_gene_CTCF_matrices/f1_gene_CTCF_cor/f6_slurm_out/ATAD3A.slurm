#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATAD3A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATAD3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATAD3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATAD3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATAD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATCAY 
