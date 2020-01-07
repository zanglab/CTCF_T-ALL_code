#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PIK3R4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3R4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3R5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3R6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIKFYVE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PILRA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PILRB 
