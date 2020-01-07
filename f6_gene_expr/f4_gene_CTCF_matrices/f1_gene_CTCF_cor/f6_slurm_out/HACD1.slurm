#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HACD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HACD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HACD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HACD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HACD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HACE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HACL1 
