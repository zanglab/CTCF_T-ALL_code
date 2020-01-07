#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR4F21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4F21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4F3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4F4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4F5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4K13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR4K17 
