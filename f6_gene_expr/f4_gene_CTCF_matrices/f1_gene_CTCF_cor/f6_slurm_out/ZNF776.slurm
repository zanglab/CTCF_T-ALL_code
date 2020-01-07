#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF776.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF776 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF777 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF778 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF780A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF780B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF781 
