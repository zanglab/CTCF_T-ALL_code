#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TPPP3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TPPP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPRA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPRG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPRG1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPRKB 
