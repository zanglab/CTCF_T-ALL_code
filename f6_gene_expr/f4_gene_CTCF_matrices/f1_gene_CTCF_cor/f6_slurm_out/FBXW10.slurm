#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FBXW10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXW10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXW11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXW12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXW2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXW4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXW4P1 
