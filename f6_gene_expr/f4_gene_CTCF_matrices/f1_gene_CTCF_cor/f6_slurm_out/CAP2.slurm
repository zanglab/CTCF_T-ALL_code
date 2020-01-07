#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CAP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPN10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPN11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPN12 
