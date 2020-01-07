#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM170A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM170A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM170B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM171 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM173 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM174 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM175 
