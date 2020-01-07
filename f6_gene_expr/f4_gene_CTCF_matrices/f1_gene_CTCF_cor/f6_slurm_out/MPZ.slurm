#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MPZ.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MPZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPZL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPZL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPZL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRAP 
