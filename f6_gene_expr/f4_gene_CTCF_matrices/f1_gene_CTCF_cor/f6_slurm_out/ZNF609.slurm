#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF609.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF609 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF610 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF611 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF613 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF614 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF615 
