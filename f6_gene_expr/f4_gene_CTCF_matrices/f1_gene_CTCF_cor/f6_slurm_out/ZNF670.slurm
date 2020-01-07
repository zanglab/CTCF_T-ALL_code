#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF670.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF670 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF670-ZNF695 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF671 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF672 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF674 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF675 
