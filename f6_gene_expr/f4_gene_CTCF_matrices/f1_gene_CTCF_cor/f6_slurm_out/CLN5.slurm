#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLN5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLN5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLN6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLN8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLNK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLNS1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLOCK 
