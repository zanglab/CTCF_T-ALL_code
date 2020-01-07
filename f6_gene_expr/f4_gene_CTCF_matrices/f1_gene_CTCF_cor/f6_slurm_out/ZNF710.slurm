#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF710.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF710 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF711 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF713 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF714 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF716 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF717 
