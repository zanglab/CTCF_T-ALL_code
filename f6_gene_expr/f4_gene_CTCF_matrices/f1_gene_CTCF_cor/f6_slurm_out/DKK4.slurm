#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DKK4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DKK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DKKL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLEC1 
