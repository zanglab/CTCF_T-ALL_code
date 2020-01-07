#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EMC3-AS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EMC3-AS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMC8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMC9 
