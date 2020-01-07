#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF558.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF558 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF559 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF559-ZNF177 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF560 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF561 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF562 
