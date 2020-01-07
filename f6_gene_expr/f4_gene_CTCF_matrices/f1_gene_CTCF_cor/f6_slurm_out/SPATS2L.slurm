#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPATS2L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPATS2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPC24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPC25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPCS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPCS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPCS3 
