#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF644.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF644 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF646 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF648 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF649 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF652 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF653 
