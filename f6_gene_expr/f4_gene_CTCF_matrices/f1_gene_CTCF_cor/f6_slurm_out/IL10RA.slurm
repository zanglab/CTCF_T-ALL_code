#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IL10RA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IL10RA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL10RB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL11RA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL12A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL12B 
