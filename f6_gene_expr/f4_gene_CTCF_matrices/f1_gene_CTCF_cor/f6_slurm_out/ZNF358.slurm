#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF358.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF358 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF362 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF365 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF366 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF367 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF37A 
