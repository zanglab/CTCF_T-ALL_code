#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C20orf24.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C20orf24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C20orf27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C20orf85 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C20orf96 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C21orf140 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C21orf2 
