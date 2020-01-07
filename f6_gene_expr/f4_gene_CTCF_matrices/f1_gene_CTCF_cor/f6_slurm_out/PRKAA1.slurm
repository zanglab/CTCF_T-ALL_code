#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRKAA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKAA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKAA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKAB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKAB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKACA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKACB 
