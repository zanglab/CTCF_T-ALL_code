#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CT45A2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CT45A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CT45A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CT45A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CT45A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CT45A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CT45A8 
