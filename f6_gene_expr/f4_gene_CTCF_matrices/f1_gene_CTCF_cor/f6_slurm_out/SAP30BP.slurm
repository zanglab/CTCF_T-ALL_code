#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SAP30BP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SAP30BP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAP30L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAPCD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAPCD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAR1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAR1B 
