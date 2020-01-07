#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR51M1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51M1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51Q1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51S1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51T1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR51V1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52A1 
