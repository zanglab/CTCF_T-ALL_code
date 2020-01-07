#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF814.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF814 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF815P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF816 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF816-ZNF321P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF821 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF823 
