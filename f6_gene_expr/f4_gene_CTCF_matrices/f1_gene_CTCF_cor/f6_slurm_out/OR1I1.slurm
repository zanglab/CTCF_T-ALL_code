#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR1I1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1I1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1J1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1J2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1J4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1K1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR1L1 
