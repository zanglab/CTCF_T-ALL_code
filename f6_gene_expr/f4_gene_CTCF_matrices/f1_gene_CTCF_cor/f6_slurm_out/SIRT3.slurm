#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SIRT3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SIRT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIRT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIRT5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIRT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIRT7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIT1 
