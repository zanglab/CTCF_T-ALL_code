#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF654.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF654 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF655 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF658 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF658B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF660 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF662 
