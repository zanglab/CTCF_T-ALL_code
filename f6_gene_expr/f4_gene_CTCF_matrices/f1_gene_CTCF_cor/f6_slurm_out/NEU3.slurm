#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NEU3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NEU3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEU4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEURL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEURL1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEURL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEURL3 
