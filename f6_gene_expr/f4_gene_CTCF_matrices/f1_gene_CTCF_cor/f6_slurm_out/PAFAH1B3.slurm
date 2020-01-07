#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PAFAH1B3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PAFAH1B3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAFAH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAGE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAGE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAGE2B 
