#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM229B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM229B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM230 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM231 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM232 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM233 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM234 
