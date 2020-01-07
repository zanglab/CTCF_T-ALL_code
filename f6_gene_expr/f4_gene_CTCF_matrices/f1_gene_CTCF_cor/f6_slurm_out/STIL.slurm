#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STIL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STIL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STIM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STIM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STK10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STK11 
