#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF663P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF663P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF664 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF665 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF667 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF668 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF669 
