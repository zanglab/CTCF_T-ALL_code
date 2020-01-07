#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PITPNM3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PITPNM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PITRM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PITX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PITX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PITX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIWIL1 
