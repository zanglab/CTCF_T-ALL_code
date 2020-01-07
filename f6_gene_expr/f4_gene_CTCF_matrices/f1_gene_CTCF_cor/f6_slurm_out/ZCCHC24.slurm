#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZCCHC24.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCCHC24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCCHC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCCHC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCCHC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCCHC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCCHC8 
