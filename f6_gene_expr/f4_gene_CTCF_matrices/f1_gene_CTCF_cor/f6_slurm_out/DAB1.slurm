#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DAB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DAB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAB2IP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DACH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DACH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DACT1 
