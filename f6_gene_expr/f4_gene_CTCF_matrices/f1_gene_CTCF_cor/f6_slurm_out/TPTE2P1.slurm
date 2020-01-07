#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TPTE2P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TPTE2P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPTE2P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPTE2P5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPTE2P6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPTEP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPX2 
