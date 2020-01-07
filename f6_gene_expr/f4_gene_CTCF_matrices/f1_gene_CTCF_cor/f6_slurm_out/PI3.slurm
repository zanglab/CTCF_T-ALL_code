#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PI3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PI3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PI4K2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PI4K2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PI4KA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PI4KAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PI4KAP2 
