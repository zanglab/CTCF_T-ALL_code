#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF239.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF239 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF248 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF250 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF251 
