#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PDRG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PDRG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDS5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDS5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDSS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDSS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDX1 
