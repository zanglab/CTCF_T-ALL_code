#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRT5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT6A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT6C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT71 
