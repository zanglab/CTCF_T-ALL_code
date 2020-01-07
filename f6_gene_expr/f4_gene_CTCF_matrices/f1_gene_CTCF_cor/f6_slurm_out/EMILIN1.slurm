#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EMILIN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EMILIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMILIN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMILIN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EML1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EML2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EML3 
