#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NACC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NACC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NACC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NADK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NADK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NADSYN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAE1 
