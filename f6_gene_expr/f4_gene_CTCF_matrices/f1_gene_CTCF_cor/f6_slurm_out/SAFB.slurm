#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SAFB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SAFB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAFB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAGE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SALL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SALL2 
