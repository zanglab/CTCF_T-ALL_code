#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VTA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VTA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VTCN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VTI1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VTI1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VWA1 
