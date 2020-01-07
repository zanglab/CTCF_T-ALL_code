#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRG3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRICKLE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRICKLE2 
