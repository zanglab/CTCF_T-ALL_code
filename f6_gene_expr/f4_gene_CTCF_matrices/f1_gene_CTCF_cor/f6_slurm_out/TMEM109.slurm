#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM109.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM109 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM110-MUSTN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM114 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM115 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM116 
