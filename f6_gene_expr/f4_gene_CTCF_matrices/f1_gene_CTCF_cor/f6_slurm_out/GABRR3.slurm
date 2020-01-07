#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GABRR3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GADD45A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GADD45B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GADD45G 
