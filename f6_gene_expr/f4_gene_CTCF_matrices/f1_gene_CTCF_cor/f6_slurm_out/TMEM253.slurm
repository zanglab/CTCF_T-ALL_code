#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM253.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM253 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM254 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM255A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM255B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM256 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM256-PLSCR3 
