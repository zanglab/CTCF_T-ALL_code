#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF318.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF318 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF319 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF320 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF321P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF322 
