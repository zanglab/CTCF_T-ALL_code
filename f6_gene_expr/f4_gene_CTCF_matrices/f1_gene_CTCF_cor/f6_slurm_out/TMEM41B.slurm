#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM41B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM41B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM45A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM45B 
