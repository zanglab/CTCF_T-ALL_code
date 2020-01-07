#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NUP153.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP153 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP155 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP160 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP188 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP205 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP210 
