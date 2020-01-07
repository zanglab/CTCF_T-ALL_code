#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM191B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM191B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM191C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM192 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM196 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM198 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM198B 
