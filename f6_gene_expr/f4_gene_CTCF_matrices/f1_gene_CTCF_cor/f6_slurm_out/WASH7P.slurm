#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WASH7P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WASH7P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WASL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WBP11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WBP11P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WBP1L 
