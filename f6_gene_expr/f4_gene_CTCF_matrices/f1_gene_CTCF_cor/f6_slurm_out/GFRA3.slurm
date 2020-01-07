#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GFRA3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GFRA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GFRA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GFRAL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GFY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGA2 
