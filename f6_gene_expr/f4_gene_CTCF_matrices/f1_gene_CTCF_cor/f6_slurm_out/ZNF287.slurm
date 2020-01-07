#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF287.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF287 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF292 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF296 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF300 
