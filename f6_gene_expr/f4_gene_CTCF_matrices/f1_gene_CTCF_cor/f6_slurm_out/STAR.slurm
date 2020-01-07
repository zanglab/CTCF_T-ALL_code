#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STAR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STAR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STARD10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STARD13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STARD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STARD3NL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STARD4 
