#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF718.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF718 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF720 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF721 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF726 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF727 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF728 
