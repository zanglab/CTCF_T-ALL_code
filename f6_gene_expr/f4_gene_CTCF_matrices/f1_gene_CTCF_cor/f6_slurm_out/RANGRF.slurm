#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RANGRF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RANGRF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAP1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAP1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAP1GAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAP1GAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAP1GDS1 
