#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR14J1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR14J1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1D4 
