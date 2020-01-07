#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF518A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF518A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF518B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF519 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF521 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF524 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF525 
