#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRT19P2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT19P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT222 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT24 
