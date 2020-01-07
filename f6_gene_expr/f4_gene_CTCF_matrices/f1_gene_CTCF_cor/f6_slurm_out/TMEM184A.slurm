#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM184A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM184A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM184B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM184C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM185A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM185B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM186 
