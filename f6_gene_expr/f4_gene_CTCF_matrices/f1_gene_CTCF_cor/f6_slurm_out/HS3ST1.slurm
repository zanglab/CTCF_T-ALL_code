#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HS3ST1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HS3ST1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HS3ST2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HS3ST3A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HS3ST3B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HS3ST4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HS3ST5 
