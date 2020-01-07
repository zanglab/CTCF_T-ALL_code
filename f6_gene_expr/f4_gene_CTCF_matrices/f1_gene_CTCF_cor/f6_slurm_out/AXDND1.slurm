#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AXDND1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AXDND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AXIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AXIN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AXL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AZGP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AZGP1P1 
