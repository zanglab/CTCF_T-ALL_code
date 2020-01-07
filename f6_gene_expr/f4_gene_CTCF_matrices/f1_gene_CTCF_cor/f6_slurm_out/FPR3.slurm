#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FPR3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FPR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRA10AC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRAS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FREM1 
