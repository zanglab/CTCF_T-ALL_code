#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o YIF1A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i YIF1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YIF1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YIPF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YIPF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YIPF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YIPF4 
