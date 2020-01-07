#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLITRK1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLITRK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLITRK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLITRK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLITRK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLITRK5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLITRK6 
