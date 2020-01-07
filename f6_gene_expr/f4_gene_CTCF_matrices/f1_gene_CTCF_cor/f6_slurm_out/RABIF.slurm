#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RABIF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RABIF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RABL2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RABL2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RABL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RABL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAC1 
