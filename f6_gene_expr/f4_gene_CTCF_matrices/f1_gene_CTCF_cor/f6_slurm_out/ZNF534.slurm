#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF534.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF534 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF536 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF540 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF541 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF542P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF543 
