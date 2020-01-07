#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPTAN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPTAN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPTB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPTBN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPTBN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPTBN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPTBN5 
