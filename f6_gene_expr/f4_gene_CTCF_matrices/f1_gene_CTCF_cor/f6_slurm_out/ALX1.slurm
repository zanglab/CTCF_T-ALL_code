#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ALX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ALX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALYREF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMACR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMBN 
