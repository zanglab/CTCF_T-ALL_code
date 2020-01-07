#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHST15.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHST15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHST2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHST3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHST4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHST5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHST6 
