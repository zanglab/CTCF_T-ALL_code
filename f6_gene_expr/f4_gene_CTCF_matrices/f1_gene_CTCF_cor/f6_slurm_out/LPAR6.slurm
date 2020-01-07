#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LPAR6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LPAR6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPCAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPCAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPCAT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPCAT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LPGAT1 
