#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DGCR6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DGCR6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGCR6L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGCR8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGKA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGKB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGKD 
