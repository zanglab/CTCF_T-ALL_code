#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM258.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM258 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM259 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM260 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM262 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM263 
