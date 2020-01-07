#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF100.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF100 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF101 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF106 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF107 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF112 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF114 
