#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o APOBEC3H.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i APOBEC3H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOBEC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOBR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOC1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOC2 
