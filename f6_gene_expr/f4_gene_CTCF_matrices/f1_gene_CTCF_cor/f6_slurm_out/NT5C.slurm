#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NT5C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NT5C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NT5C1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NT5C1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NT5C1B-RDH14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NT5C2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NT5C3A 
