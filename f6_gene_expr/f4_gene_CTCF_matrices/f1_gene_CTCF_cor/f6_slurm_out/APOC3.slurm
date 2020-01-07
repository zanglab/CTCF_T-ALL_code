#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o APOC3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i APOC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOC4-APOC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOF 
