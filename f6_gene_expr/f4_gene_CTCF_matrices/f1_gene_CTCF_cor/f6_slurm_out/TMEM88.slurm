#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM88.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM88 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM88B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM89 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM8A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM8B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM9 
