#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BECN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BECN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BECN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BEGAIN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BEND2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BEND3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BEND3P3 
