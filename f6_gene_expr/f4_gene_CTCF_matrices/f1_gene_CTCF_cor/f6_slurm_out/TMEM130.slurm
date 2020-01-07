#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM130.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM130 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM131 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM132A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM132B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM132C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM132D 
