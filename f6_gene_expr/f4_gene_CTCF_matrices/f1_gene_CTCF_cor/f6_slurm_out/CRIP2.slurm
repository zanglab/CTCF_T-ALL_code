#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CRIP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CRIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRIP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRIPT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRISP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRISP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRISP3 
