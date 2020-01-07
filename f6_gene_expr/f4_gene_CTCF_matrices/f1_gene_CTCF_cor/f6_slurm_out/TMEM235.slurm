#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM235.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM235 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM236 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM237 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM238 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM239 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM240 
