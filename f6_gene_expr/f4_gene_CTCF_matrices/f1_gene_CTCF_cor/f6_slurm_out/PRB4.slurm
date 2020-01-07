#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRB4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRCC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRCD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRCP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDM1 
