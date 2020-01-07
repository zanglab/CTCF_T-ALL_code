#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM140.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM140 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM141 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM143 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM144 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM145 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM147 
