#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EDAR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EDAR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDARADD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDDM3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDDM3B 
