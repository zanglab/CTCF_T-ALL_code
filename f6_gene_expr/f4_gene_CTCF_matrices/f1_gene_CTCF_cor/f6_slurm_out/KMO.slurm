#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KMO.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KMO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KMT2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KMT2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KMT2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KMT2D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KMT2E 
