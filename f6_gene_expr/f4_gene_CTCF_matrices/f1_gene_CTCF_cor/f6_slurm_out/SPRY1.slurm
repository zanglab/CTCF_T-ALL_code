#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPRY1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPRY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPRY2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPRY3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPRY4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPRYD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPRYD4 
