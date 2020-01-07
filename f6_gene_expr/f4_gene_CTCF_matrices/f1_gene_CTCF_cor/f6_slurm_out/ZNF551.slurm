#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF551.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF551 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF552 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF554 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF555 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF556 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF557 
