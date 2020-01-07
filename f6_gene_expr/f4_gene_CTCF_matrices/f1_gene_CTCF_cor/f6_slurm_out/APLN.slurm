#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o APLN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i APLN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APLNR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APLP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APLP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APMAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOA1 
