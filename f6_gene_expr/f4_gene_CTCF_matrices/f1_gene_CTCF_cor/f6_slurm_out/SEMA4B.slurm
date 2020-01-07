#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SEMA4B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SEMA4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEMA4C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEMA4D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEMA4F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEMA4G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEMA5A 
