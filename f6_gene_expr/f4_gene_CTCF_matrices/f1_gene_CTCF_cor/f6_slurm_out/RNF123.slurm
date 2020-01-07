#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNF123.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF123 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF125 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF126 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF126P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF128 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF13 
