#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TSPAN9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPEAR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPY10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPY2 
