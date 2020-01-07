#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZCCHC12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCCHC12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCCHC13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCCHC14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCCHC17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCCHC18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCCHC2 
