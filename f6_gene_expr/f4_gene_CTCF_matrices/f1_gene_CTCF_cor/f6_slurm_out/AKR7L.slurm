#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AKR7L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AKR7L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKT1S1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AKTIP 
