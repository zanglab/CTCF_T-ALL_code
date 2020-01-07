#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF260.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF260 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF263 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF264 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF266 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF267 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF268 
