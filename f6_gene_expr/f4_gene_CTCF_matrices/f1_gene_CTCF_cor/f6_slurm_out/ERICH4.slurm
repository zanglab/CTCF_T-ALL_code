#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ERICH4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ERICH4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERICH5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERICH6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERICH6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERLEC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERLIN1 
