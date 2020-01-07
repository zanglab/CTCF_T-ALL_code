#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PDLIM7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PDLIM7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDPK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDPN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDPR 
