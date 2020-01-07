#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM151A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM151A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM151B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM154 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM155 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM156 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM158 
