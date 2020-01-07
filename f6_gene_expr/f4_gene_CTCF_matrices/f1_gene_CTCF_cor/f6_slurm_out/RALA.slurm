#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RALA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RALA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RALB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RALBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RALGAPA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RALGAPA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RALGAPB 
