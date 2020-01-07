#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM47.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM50A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM50B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM51 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM52B 
