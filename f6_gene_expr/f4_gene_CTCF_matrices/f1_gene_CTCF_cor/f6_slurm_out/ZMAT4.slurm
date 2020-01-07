#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZMAT4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMAT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMAT5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMIZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMIZ2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMPSTE24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMYM1 
