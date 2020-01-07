#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM91.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM91 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM92 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM95 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM97 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM98 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM99 
