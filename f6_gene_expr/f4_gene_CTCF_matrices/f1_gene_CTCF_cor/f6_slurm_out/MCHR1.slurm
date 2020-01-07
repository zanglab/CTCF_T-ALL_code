#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MCHR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MCHR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCHR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCIDAS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCM10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCM2 
