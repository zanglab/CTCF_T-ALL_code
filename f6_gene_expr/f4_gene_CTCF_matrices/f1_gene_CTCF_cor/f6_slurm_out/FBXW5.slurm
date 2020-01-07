#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FBXW5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXW5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXW7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXW8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXW9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCAMR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCAR 
