#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF888.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF888 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF890P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF891 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF90 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF91 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF92 
