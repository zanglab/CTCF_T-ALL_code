#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRT32.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT33A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT33B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT36 
