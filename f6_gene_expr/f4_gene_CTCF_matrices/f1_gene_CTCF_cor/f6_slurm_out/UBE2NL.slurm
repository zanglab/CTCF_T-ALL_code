#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UBE2NL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2NL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2O 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2Q1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2Q2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2Q2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2Q2P1 
