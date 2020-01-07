#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR2A2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2A20P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2A25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2A42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2A5 
