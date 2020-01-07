#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRT83.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT83 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT84 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT85 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT86 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT8P41 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT9 
