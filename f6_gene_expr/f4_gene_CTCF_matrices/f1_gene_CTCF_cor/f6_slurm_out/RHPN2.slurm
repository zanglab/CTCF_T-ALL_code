#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RHPN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RHPN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIBC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIBC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIC8A 
