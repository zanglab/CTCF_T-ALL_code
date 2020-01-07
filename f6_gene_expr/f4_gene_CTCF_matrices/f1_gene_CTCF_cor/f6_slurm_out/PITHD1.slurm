#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PITHD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PITHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PITPNA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PITPNB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PITPNC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PITPNM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PITPNM2 
