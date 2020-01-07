#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLURP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLURP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLX1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLX1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLX4IP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMAD1 
