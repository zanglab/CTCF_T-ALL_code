#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF37BP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF37BP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF382 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF383 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF384 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF385A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF385B 
