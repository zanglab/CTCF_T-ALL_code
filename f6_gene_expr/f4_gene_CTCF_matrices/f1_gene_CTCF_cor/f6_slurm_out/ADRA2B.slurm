#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADRA2B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADRA2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADRA2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADRB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADRB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADRB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADRM1 
