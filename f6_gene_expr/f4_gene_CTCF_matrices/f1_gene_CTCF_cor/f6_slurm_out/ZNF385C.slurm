#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF385C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF385C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF385D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF391 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF394 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF395 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF396 
