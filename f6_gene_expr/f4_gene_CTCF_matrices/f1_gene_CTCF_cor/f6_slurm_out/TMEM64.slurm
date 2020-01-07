#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM64.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM64 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM65 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM67 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM68 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM69 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM70 
