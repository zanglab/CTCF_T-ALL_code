#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF865.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF865 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF876P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF878 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF879 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF880 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF883 
