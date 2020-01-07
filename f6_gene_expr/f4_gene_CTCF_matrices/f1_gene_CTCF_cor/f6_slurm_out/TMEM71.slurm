#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM71.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM71 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM72 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM74 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM74B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM79 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM80 
