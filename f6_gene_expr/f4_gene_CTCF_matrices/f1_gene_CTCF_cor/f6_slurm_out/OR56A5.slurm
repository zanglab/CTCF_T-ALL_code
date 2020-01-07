#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR56A5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR56A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR56B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR56B4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR5AK2 
