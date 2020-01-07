#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRIM50.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM50 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM51 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM53AP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM55 
