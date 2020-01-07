#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RIOK2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RIOK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIOK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIPK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIPK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIPK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIPK4 
