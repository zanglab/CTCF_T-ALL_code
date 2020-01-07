#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IL20.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IL20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL20RA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL20RB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL21R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL22 
