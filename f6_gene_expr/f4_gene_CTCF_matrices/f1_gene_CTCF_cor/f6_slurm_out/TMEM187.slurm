#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM187.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM187 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM189 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM189-UBE2V1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM190 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM191A 
