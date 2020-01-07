#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UBE2G1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2G1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2G2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2I 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2J1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2J2 
