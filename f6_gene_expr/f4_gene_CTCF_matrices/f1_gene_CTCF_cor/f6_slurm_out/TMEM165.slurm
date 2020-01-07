#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM165.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM165 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM167A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM167B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM168 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM169 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM17 
