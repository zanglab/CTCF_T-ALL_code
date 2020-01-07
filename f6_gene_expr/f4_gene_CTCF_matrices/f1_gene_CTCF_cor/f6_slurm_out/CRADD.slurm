#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CRADD.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CRADD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRBN 
