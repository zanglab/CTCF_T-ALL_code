#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF593.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF593 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF594 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF595 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF596 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF597 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF598 
