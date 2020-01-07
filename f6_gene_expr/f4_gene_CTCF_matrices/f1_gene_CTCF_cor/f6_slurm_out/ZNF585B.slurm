#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF585B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF585B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF586 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF587 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF587B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF589 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF592 
