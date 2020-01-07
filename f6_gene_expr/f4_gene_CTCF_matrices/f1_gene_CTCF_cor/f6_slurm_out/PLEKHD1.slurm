#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLEKHD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHG3 
