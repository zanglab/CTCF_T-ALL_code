#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF410.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF410 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF414 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF415 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF416 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF417 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF418 
