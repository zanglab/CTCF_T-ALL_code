#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF117.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF117 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF121 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF124 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF131 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF132 
