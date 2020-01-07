#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VWA2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VWA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VWA3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VWA3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VWA5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VWA5B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VWA5B2 
