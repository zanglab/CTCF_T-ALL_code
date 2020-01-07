#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UBTFL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UBTFL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBXN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBXN10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBXN11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBXN2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBXN2B 
