#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATRAID.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATRAID 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATRIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATRN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATRNL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATRX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATXN1 
