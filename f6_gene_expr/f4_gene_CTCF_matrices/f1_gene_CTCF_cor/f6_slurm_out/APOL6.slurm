#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o APOL6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i APOL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOLD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOOL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOOP5 
