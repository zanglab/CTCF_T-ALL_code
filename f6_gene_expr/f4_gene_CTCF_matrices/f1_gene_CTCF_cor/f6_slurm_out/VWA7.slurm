#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VWA7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VWA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VWA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VWC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VWC2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VWCE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VWDE 
