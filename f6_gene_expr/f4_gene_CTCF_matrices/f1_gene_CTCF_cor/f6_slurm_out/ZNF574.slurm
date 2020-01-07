#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF574.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF574 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF575 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF576 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF577 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF578 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF579 
