#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TTLL7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TTLL7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTLL9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTPAL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTR 
