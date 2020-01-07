#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SIN3B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SIN3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIPA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIPA1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIPA1L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIPA1L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIRPA 
