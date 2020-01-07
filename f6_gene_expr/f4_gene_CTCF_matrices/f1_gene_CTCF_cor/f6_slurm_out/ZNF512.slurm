#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF512.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF512 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF512B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF513 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF514 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF516 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF517 
