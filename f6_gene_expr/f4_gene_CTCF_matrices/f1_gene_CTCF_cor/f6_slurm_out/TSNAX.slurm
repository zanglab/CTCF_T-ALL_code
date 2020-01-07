#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TSNAX.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TSNAX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSNAX-DISC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSNAXIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN11 
