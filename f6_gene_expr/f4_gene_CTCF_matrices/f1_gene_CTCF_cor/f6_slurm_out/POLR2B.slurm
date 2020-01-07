#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POLR2B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2G 
