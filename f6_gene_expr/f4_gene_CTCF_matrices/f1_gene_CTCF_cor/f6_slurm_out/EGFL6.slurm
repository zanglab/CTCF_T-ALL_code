#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EGFL6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EGFL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EGFL7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EGFL8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EGFLAM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EGFR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EGLN1 
