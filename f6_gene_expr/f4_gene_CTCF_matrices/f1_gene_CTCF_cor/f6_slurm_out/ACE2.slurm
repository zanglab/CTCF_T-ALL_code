#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACE2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACER1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACER2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACER3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACHE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACIN1 
