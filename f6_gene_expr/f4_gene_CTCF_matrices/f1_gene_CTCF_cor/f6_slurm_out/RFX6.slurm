#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RFX6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RFX6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFX7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFX8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFXANK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFXAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGCC 
