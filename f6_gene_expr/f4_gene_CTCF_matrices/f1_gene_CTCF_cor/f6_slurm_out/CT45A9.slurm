#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CT45A9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CT45A9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CT47A12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CT47B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CT55 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CT62 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CT83 
