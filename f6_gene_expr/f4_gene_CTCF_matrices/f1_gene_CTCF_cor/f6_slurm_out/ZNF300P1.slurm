#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF300P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF300P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF302 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF304 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF311 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF316 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF317 
