#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PADI4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PADI4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PADI6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAEP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAFAH1B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAFAH1B2 
