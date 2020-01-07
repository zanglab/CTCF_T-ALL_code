#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACO1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACOT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACOT11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACOT12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACOT13 
