#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CMAS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CMAS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMBL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMIP 
