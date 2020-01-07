#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RFWD3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RFWD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFX5 
