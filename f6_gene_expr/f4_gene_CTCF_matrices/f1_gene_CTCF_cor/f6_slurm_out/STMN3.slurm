#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STMN3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STMN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STMN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STMND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STOM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STOML1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STOML2 
