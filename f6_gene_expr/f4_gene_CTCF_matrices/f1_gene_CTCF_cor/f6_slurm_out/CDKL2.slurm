#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDKL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDKL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDKL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDKL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDKL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDKN1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDKN1B 
