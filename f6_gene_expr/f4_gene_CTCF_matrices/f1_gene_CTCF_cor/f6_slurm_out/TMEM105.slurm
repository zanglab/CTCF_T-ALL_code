#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM105.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM105 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM106A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM106B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM106C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM107 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM108 
