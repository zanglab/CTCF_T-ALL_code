#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DGKZ.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DGKZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGUOK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHCR24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHCR7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHDDS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHDH 
