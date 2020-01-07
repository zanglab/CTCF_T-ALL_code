#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SYNGR2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNGR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNGR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNGR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNJ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNJ2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNJ2BP 
