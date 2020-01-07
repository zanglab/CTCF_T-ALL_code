#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ERI1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ERI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERI3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERICH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERICH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERICH3 
