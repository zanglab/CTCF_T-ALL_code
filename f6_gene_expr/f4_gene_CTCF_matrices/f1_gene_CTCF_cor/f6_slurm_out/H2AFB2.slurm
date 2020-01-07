#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o H2AFB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i H2AFB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H2AFB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H2AFJ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H2AFV 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H2AFX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H2AFY 
