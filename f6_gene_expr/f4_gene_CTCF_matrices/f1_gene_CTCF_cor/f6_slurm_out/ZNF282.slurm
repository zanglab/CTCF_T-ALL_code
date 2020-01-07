#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF282.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF282 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF283 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF284 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF285 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF286A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF286B 
