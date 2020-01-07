#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMCO3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMCO3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMCO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMCO5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMCO5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMCO6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMED1 
