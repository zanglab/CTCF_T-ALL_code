#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPAG4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPAG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPAG5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPAG6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPAG7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPAG8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPAG9 
