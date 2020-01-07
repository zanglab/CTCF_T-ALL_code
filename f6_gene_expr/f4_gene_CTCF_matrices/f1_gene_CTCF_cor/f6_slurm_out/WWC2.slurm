#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WWC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WWC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WWC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WWOX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WWP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WWP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WWTR1 
