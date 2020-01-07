#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IL22RA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IL22RA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL22RA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL23A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL23R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL25 
