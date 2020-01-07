#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PIK3CD.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3CD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3CG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3IP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3R1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3R2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3R3 
