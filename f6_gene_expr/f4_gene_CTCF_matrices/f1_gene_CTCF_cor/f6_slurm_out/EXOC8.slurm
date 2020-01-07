#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EXOC8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOC8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOSC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOSC10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOSC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOSC3 
