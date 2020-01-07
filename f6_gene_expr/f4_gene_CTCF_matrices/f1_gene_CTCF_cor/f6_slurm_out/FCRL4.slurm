#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FCRL4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FCRL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCRL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCRL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCRLA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCRLB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FDCSP 
