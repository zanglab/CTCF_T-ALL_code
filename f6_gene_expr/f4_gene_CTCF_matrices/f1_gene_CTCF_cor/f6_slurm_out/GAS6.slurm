#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GAS6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GAS6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAS7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAS8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAST 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GATA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GATA2 
