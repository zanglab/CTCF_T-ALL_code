#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM132E.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM132E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM134 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM135 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM136 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM138 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM139 
