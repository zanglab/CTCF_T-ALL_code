#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPEM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPEM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPEN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPERT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPESP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPG11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPG21 
