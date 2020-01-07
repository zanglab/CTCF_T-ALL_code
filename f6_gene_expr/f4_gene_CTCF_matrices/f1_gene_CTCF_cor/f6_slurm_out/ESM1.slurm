#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ESM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ESM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ESPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ESPN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ESPNL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ESPNP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ESR1 
