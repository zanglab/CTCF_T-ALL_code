#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF436.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF436 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF438 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF439 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF440 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF441 
