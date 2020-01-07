#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AKR1C8P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR1C8P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR1D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR1E2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR7A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR7A2P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR7A3 
