#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NABP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NABP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NABP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NACA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NACA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NACAD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NACAP1 
