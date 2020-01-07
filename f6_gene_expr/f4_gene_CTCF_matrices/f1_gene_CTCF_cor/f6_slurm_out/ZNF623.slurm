#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF623.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF623 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF624 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF625 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF625-ZNF20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF626 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF627 
