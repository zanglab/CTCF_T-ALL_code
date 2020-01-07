#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o USO1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i USO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP13 
