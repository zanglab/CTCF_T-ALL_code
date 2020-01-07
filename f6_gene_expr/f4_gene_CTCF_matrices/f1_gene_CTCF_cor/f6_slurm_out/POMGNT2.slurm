#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POMGNT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POMGNT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POMK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POMP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POMT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POMT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POMZP3 
