#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF324.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF324 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF324B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF326 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF329 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF330 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF331 
