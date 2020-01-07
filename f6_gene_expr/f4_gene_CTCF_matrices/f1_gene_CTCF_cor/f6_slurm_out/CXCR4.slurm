#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CXCR4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CXCR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXCR5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXCR6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXXC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXXC1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXXC4 
