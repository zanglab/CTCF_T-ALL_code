#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o E2F3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i E2F3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i E2F4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i E2F5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i E2F6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i E2F7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i E2F8 
