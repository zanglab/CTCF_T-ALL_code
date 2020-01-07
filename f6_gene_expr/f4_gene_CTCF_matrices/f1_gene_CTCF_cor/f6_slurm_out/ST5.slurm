#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ST5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ST5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST6GAL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST6GAL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST6GALNAC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST6GALNAC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST6GALNAC3 
