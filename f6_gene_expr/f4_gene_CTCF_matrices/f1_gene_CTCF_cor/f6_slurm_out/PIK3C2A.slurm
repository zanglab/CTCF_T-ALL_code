#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PIK3C2A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3C2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3C2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3C2G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3C3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3CA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIK3CB 
