#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UBE2B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2D3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2D4 
