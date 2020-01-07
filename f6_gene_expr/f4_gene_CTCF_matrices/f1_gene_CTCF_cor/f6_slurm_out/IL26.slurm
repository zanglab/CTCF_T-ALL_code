#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IL26.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IL26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL27RA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL2RA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL2RB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL2RG 
