#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NDC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDC80 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDEL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDFIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDFIP2 
