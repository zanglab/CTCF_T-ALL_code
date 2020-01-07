#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZSWIM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSWIM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSWIM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSWIM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSWIM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSWIM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSWIM6 
