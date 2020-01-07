#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF682.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF682 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF683 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF684 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF687 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF688 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF689 
