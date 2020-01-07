#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SERPINI2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERTAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERTAD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERTAD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERTAD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERTM1 
