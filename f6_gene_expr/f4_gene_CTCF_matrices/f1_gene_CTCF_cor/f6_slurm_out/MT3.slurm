#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MT3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTAP 
