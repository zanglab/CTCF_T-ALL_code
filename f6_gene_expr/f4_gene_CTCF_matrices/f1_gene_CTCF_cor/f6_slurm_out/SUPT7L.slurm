#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SUPT7L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SUPT7L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUPV3L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SURF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SURF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SURF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SURF6 
