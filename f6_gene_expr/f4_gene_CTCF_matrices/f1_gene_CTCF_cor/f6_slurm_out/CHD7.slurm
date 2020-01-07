#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHD7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHD9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHDH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHEK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHEK2 
