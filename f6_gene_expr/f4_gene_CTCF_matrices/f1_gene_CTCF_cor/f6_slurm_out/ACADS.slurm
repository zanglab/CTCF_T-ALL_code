#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACADS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACADS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACADSB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACADVL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACAN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACAP2 
