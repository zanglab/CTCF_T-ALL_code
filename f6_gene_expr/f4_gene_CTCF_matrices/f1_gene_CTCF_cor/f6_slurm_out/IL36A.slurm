#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IL36A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IL36A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL36B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL36G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL36RN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL3RA 
