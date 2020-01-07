#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF544.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF544 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF546 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF547 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF548 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF549 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF550 
