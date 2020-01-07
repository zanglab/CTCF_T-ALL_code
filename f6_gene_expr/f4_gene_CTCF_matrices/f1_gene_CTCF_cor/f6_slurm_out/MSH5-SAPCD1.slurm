#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MSH5-SAPCD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MSH5-SAPCD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSH6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSL2 
