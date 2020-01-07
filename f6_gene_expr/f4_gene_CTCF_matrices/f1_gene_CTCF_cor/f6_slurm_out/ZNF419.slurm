#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF419.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF419 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF420 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF423 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF425 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF426 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF428 
