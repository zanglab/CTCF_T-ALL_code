#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o A1BG.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i A1BG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i A1CF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i A2M 
time python ../7_cor_CTCF_binding_GeneExpr.py -i A2ML1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i A2MP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i A3GALT2 
