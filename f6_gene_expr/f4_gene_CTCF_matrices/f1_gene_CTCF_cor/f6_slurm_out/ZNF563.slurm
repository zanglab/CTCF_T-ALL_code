#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF563.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF563 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF564 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF565 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF566 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF567 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF568 
