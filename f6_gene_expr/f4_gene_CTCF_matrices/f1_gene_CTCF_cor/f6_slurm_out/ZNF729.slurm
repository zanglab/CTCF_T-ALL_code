#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF729.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF729 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF730 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF732 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF735 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF736 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF737 
