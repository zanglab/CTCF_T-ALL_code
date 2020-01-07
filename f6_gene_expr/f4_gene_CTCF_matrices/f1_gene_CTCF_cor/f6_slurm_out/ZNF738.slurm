#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF738.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF738 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF74 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF740 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF746 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF747 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF749 
