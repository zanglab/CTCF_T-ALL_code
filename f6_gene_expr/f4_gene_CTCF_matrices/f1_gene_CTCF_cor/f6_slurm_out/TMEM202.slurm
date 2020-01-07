#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM202.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM202 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM203 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM204 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM205 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM206 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM207 
