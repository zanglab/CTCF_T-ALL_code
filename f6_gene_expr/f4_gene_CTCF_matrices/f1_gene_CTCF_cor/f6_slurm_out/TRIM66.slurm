#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRIM66.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM66 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM67 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM68 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM69 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM71 
