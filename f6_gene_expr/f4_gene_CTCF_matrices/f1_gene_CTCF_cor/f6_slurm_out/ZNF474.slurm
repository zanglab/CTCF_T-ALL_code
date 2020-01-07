#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF474.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF474 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF479 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF48 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF480 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF483 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF484 
